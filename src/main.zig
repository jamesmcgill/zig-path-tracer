const c_import = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");
const math = @import("std").math;
const time = @import("std").time;

const Vec3 = @import("Vec3.zig").Vec3;
const Ray = @import("Ray.zig").Ray;

//------------------------------------------------------------------------------
const Vec2 = struct
{
  x: f32,
  y: f32,
};

//------------------------------------------------------------------------------
// Return a random number in the range [0, 1)
//------------------------------------------------------------------------------
const MyRand = struct
{
  rng: std.rand.DefaultPrng,

  //----------------------------------------------------------------------------
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

  //----------------------------------------------------------------------------
  pub fn float(self: *MyRand) f32
  {
    return self.rng.random.float(f32);
  }

  //----------------------------------------------------------------------------
  pub fn randomPointFromUnitSphere(self: *MyRand) Vec3
  {
    return while (true)
    {
      const p = Vec3 {
        .x = (self.float() * 2.0) - 1.0,
        .y = (self.float() * 2.0) - 1.0,
        .z = (self.float() * 2.0) - 1.0,
      };
      if (p.lengthSq() < 1.0) { break p; }
    } else p;
  }
  //----------------------------------------------------------------------------

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
const Camera = struct
{
  position: Vec3,
  scale: Vec2,
  offset: Vec2,
  frustum_dist: f32,

  const default_frustum_width: f32 = 200.0;
  const default_frustum_height: f32 = 200.0;
  const default_frustum_dist: f32 = 100.0;

  //----------------------------------------------------------------------------
  // Scale image (screen) to world space
  //----------------------------------------------------------------------------
  pub fn calcScale(
    image_width: u32,
    image_height: u32,
    frustum_width: f32,
    frustum_height: f32) Vec2
  {
    return Vec2
    {
      .x = frustum_width / @intToFloat(f32, image_width - 1),
      .y = frustum_height / @intToFloat(f32, image_height - 1),
    };
  }

  //----------------------------------------------------------------------------
  // Center image at x=0, y=0
  //----------------------------------------------------------------------------
  pub fn calcOffset(frustum_width: f32, frustum_height: f32) Vec2
  {
    return Vec2
    {
      .x = -frustum_width / 2.0,
      .y = -frustum_height / 2.0,
    };
  }

  //----------------------------------------------------------------------------
  pub fn create(
    position: Vec3,
    image_width: u32,
    image_height: u32,
    frustum_width: f32,
    frustum_height: f32,
    frustum_dist: f32) Camera
  {
    return Camera
    {
        .position = position,
        .scale = calcScale(
          image_width, image_height, frustum_width, frustum_height),
        .offset = calcOffset(frustum_width, frustum_height),
        .frustum_dist = frustum_dist,
    };
  }

  //----------------------------------------------------------------------------
  pub fn calcRay(cam: Camera, col: f32, row: f32, rand: *MyRand) Ray
  {
    // The point that corresponds to on the frustum plane
    const to_x: f32 = cam.offset.x + ((col + rand.float()) * cam.scale.x);
    const to_y: f32 = cam.offset.y + ((row + rand.float()) * cam.scale.y);
    const to_pixel = Vec3{.x = to_x, .y = to_y, .z = cam.frustum_dist};

    const ray_dir: Vec3 = to_pixel.subtract(cam.position);
    return Ray{.origin = cam.position, .direction = ray_dir};
  }

  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
const HitInfo = struct
{
  point: Vec3,
  t: f32,
  surface_normal: Vec3,
  material_ptr: *const Material,
};

//------------------------------------------------------------------------------
const Scene = struct
{
  spheres: []const Sphere,
};

//------------------------------------------------------------------------------
const basic_scene = Scene
{
  .spheres = &[_]Sphere
  {
    .{
      .position = .{.x = 0.0, .y = -320.0, .z = 110.0},
      .radius = 300.0,
      .material =
      .{
        .Lambertian = .{.albedo = .{.x = 0.8, .y = 0.8, .z = 0.0}}
      }
    },
    .{
      .position = .{.x = 21.0, .y = 0.0, .z = 110.0},
      .radius = 20.0,
      .material =
      .{
        .Lambertian = .{.albedo = .{.x = 0.8, .y = 0.3, .z = 0.3}}
      }
    },
    .{
      .position = .{.x = -21.0, .y = 0.0, .z = 110.0},
      .radius = 20.0,
      .material =
      .{
        .Metal = .{.albedo = .{.x = 0.8, .y = 0.6, .z = 0.2}}
      }
    },
  },
};

//------------------------------------------------------------------------------
const Sphere = struct
{
  position: Vec3,
  radius: f32,
  material: Material,

  //----------------------------------------------------------------------------
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
        .material_ptr = &self.material,
      };
    }
    return null;
  }

  //----------------------------------------------------------------------------
  pub fn isInRange(val: f32, min: f32, max: f32) bool
  {
    if (val < min) { return false; }
    if (val > max) { return false; }
    return true;
  }

  //----------------------------------------------------------------------------
  pub fn closestValidT(t1: f32, t2: f32, t_min: f32, t_max: f32) ?f32
  {
    var isT1Valid = isInRange(t1, t_min, t_max);
    var isT2Valid = isInRange(t2, t_min, t_max);

    if (!isT1Valid and !isT2Valid) { return null; }
    if (isT1Valid and !isT2Valid) { return t1; }
    if (isT2Valid and !isT1Valid) { return t2; }

    return if (t1 < t2) t1 else t2;
  }

  //----------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
const MaterialTag = enum
{
  Lambertian,
  Metal,
};

const Material = union(MaterialTag)
{
  Lambertian: LambertianMaterial,
  Metal: MetalMaterial,
};

//------------------------------------------------------------------------------
const ScatterInfo = struct
{
  ray: Ray,
  attenuation: Vec3,
};

//------------------------------------------------------------------------------
const LambertianMaterial = struct
{
  albedo: Vec3,

  pub fn scatter(
    self: LambertianMaterial,
    ray: Ray,
    hit: HitInfo,
    rand: *MyRand,
  ) ?ScatterInfo
  {
    const normal_sphere_pos =  hit.point.add(hit.surface_normal.normalized());
    const reflect_point =
      normal_sphere_pos.add(rand.randomPointFromUnitSphere());

    return ScatterInfo
    {
      .ray = Ray
      {
        .origin = hit.point,
        .direction = reflect_point.subtract(hit.point),
      },
      .attenuation = self.albedo,
    };

  }
};

//------------------------------------------------------------------------------
const MetalMaterial = struct
{
  albedo: Vec3,

  pub fn scatter(
    self: MetalMaterial,
    ray: Ray,
    hit: HitInfo,
    rand: *MyRand,
  ) ?ScatterInfo
  {
    const reflected_dir = ray.direction.reflect(hit.surface_normal);
    if (reflected_dir.dot(hit.surface_normal) > 0.0)
    {
      return ScatterInfo
      {
        .ray = Ray{.origin = hit.point, .direction = reflected_dir},
        .attenuation = self.albedo,
      };
    }

    return null;
  }
};

//------------------------------------------------------------------------------
pub fn calcColor(ray: Ray, scene: *const Scene, rand: *MyRand, call_depth: u32) Vec3
{
  if (call_depth > 5) { return Vec3{.x = 0.0, .y = 1.0, .z = 1.0}; }

  // Test all objects in the scene
  var hit_something: bool = false;
  var closest_t: f32 = math.f32_max;
  var hit: HitInfo = undefined;

  for (scene.spheres) |*sphere_ptr|
  {
    if (sphere_ptr.hitTest(ray, 0.001, closest_t)) |info|
    {
      hit_something = true;
      closest_t = info.t;
      hit = info;
    }
  }

  if (hit_something)
  {
    const scatter_info = switch(hit.material_ptr.*)
    {
      MaterialTag.Lambertian => |mat| mat.scatter(ray, hit, rand),
      MaterialTag.Metal =>      |mat| mat.scatter(ray, hit, rand),
    };

    if (scatter_info) |scatter| {
      return calcColor(
        scatter.ray,
        scene,
        rand,
        call_depth + 1
      ).multiply(scatter.attenuation);
    }
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
  // Scene to draw
  const scene_ptr = &basic_scene;

  // Output image details
  const image_filename = "test.bmp";
  const image_width: u32 = 256;
  const image_height: u32 = 256;

  // Rendering parameters
  const num_samples: u32 = 100;
  const num_samples_recip: f32 = 1.0 / @intToFloat(f32, num_samples);

  var rand = MyRand.create();

  const camera_pos = Vec3{.x = 0.0, .y = 0.0, .z = 0.0};
  const camera = Camera.create(
    camera_pos,
    image_width, image_height,
    Camera.default_frustum_width,
    Camera.default_frustum_height,
    Camera.default_frustum_dist
  );

  // Output image
  var pixels: [image_width * image_height]Color = undefined;

  var timer = try time.Timer.start();
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
      const ray = camera.calcRay(col, row, &rand);
      color = color.add( calcColor(ray, scene_ptr, &rand, 0) );
    }

    color = color.scale(num_samples_recip);   // Average of samples
    color = color.sqrt();                     // Gamma correction
    item.red = @floatToInt(u8, color.x * 255.0);
    item.green = @floatToInt(u8, color.y * 255.0);
    item.blue = @floatToInt(u8, color.z * 255.0);
    item.alpha = 0xFF;
  }
  std.debug.warn("Render Time: {d:.3}s\n", .{@intToFloat(f32, timer.read()) / time.ns_per_s});

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
