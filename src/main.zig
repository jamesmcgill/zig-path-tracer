const c_import = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");
const math = @import("std").math;

const Vec3 = @import("Vec3.zig").Vec3;
const Ray = @import("Ray.zig").Ray;

//------------------------------------------------------------------------------
// Return a random number in the range [0, 1)
//------------------------------------------------------------------------------
const MyRand = struct
{
  rng: std.rand.DefaultPrng,

  pub fn create() MyRand
  {
    var buf: [8]u8 = undefined;
    std.crypto.randomBytes(buf[0..]) catch |err|
    {
      std.debug.warn("ERROR couldn't seed RNG\n", .{});
    };
    const seed = std.mem.readIntNative(u64, buf[0..8]);

    return MyRand {
      .rng = std.rand.DefaultPrng.init(seed),
    };
  }

  pub fn float(self: *MyRand) f32
  {
    return self.rng.random.float(f32);
  }
};

//------------------------------------------------------------------------------
const Color = packed struct
{
  red: u8,
  green: u8,
  blue: u8,
  alpha: u8,

  pub const NUM_COMPONENTS: u32 = 4;
};

//------------------------------------------------------------------------------
const HitInfo = struct
{
  point: Vec3,
  t: f32,
  surface_normal: Vec3,
};

//------------------------------------------------------------------------------
const Scene = struct
{
  spheres: []Sphere,
};

//------------------------------------------------------------------------------
const Sphere = struct
{
  position: Vec3,
  radius: f32,
  color: Vec3,

  pub fn hitTest(self: Sphere, ray: Ray, t_min: f32, t_max: f32) ?HitInfo
  {
    const displacement = ray.origin.subtract(self.position);

    const a: f32 = Vec3.dot(ray.direction, ray.direction);
    const b: f32 = 2.0 * Vec3.dot(displacement, ray.direction);
    const c: f32 = Vec3.dot(displacement, displacement) -
      (self.radius * self.radius);

    const discriminant: f32 = (b * b) - (4.0 * a * c);
    if (discriminant < 0.0) { return null; }

    const discrim_root = math.sqrt(discriminant);
    const double_a = 2.0 * a;
    const t1 = (-b + discrim_root) / double_a;
    const t2 = (-b - discrim_root) / double_a;

    if (closestValidT(t1, t2, t_min, t_max)) |t|
    {
      const point_at_t = ray.pointAtT(t);
      return HitInfo
      {
        .point = point_at_t,
        .t = t,
        .surface_normal = point_at_t.subtract(self.position),
      };
    }
    return null;
  }

  pub fn isInRange(val: f32, min: f32, max: f32) bool
  {
    if (val < min) { return false; }
    if (val > max) { return false; }
    return true;
  }

  pub fn closestValidT(t1: f32, t2: f32, t_min: f32, t_max: f32) ?f32
  {
    var isT1Valid = isInRange(t1, t_min, t_max);
    var isT2Valid = isInRange(t2, t_min, t_max);

    if (!isT1Valid and !isT2Valid) { return null; }
    if (isT1Valid and !isT2Valid) { return t1; }
    if (isT2Valid and !isT1Valid) { return t2; }

    return if (t1 < t2) t1 else t2;
  }
};

//------------------------------------------------------------------------------
pub fn calcColor(ray: Ray, scene: Scene) Vec3
{
  // Test all objects in the scene
  var hit_something: bool = false;
  var closest_t: f32 = 100000.0;
  var closest_hit: HitInfo = undefined;

  for (scene.spheres) |sphere|
  {
    if (sphere.hitTest(ray, 0.0, closest_t)) |info|
    {
      hit_something = true;
      closest_t = info.t;
      closest_hit = info;
    }
  }
  if (hit_something)
  {
    return closest_hit.surface_normal.normalized().rescaled();
  }

  // Draw background
  const unit_direction = ray.direction.normalized();
  const t: f32 = (unit_direction.y + 1.0) * 0.5;

  const white_tone = Vec3.init(1.0, 1.0, 1.0).scale(t);
  const blue_tone = Vec3.init(0.5, 0.7, 1.0).scale(1.0 - t);
  return white_tone.add(blue_tone);
}

//------------------------------------------------------------------------------
pub fn main() anyerror!void
{
  // Output image details
  const image_filename = "test.bmp";
  const image_width: u32 = 256;
  const image_height: u32 = 256;

  // Rendering parameters
  const num_samples: u32 = 100;
  const num_samples_recip: f32 = 1.0 / @intToFloat(f32, num_samples);

  // Scene objects
  const scene = Scene
  {
    .spheres = &[_]Sphere
    {
      .{
        .position = Vec3.init(-30.0, 0.0, 110.0),
        .radius = 20.0,
        .color = Vec3.init(1.0, 0.0, 0.0),
      },
      .{
        .position = Vec3.init(30.0, 0.0, 110.0),
        .radius = 20.0,
        .color = Vec3.init(1.0, 0.0, 0.0),
      },
    },
  };

  var r = MyRand.create();

  // Output image
  var pixels: [image_width * image_height]Color = undefined;

  // Frustum extents
  const frustum_width: f32 = 200.0;
  const frustum_height: f32 = 200.0;
  const frustum_dist: f32 = 100.0;

  // Scale image (screen) to world space
  const x_scale: f32 = frustum_width / @intToFloat(f32, image_width - 1);
  const y_scale: f32 = frustum_height / @intToFloat(f32, image_height - 1);

  // Center image at x=0, y=0
  const x_offset: f32 = -frustum_width / 2.0;
  const y_offset: f32 = -frustum_height / 2.0;

  // Ray casting from origin
  const ray_origin = Vec3 {.x = 0.0, .y = 0.0, .z = 0.0};

  for (pixels) |*item, it|
  {
    const pix = @intCast(u32, it);

    // The current pixel's coordinate position (row, col)
    const col: f32 = @intToFloat(f32, pix % image_width);
    const row: f32 = @intToFloat(f32, pix / image_width);

    var color = Vec3 {.x = 0.0, .y = 0.0, .z = 0.0};
    var samp: u32 = 0;
    while (samp < num_samples) : (samp += 1)
    {
      // The point that corresponds to on the frustum plane
      const to_x: f32 = x_offset + ((col + r.float()) * x_scale);
      const to_y: f32 = y_offset + ((row + r.float()) * y_scale);
      const to_pixel: Vec3 = Vec3.init(to_x, to_y, frustum_dist);

      const ray_dir: Vec3 = to_pixel.subtract(ray_origin);
      const ray = Ray.init(ray_origin, ray_dir);

      color = color.add( calcColor(ray, scene) );
    }

    color = color.scale(num_samples_recip);
    item.red = @floatToInt(u8, color.x * 255.0);
    item.green = @floatToInt(u8, color.y * 255.0);
    item.blue = @floatToInt(u8, color.z * 255.0);
    item.alpha = 0xFF;
  }

  // Save the image to a file
  var ret = c_import.stbi_write_bmp(
    image_filename,
    @as(c_int, image_width),
    @as(c_int, image_height),
    @as(c_int, Color.NUM_COMPONENTS),
    @ptrCast(*const c_void, &pixels[0])
  );

  if (ret == 0) {
    std.debug.warn("FAILED writing to file: {}\n", .{image_filename});
  } else {
    std.debug.warn("Successfuly written file: {}\n", .{image_filename});
  }
}

//------------------------------------------------------------------------------
