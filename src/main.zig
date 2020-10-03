const c_import = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");
const math = @import("std").math;
const time = @import("std").time;

const Vec3 = @import("Vec3.zig").Vec3;
const Ray = @import("Ray.zig").Ray;

//------------------------------------------------------------------------------
pub fn degreesToRadians(degrees: f32) f32
{
  return (degrees * math.pi) / 180.0;
}

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
  top_left: Vec3,
  horizontal: Vec3,
  vertical: Vec3,
  px_scale: Vec2,

  //----------------------------------------------------------------------------
  // Scale image (screen) to world space
  //----------------------------------------------------------------------------
  pub fn calcPxToViewport(
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
  pub fn create(
    lookfrom: Vec3,
    lookat: Vec3,
    vup: Vec3,
    image_width: u32,
    image_height: u32,
    vertical_fov: f32,
    aspect_ratio: f32) Camera
  {
    const theta = degreesToRadians(vertical_fov);
    const h = math.tan(theta / 2.0);
    const frustum_height = 2.0 * h;
    const frustum_width = aspect_ratio * frustum_height;

    const look_dir = lookfrom.subtract(lookat); // -z is forward
    const w = look_dir.normalized();
    const u = vup.cross(w).normalized();
    const v = w.cross(u);

    const half_horiz = u.scale(frustum_width / 2.0);
    const half_vert = v.scale(frustum_height / 2.0);
    const top_left_pos = lookfrom.subtract(half_horiz).add(half_vert).subtract(w);

    const px_proj = calcPxToViewport(image_width, image_height, frustum_width, frustum_height);

    return Camera
    {
      .position = lookfrom,
      .top_left = top_left_pos,
      .horizontal = u,
      .vertical = v,
      .px_scale = px_proj,
    };
  }

  //----------------------------------------------------------------------------
  pub fn calcRay(cam: Camera, col: f32, row: f32, rand: *MyRand) Ray
  {
    // The point that corresponds to on the frustum plane
    const u = cam.horizontal.scale((col + rand.float()) * cam.px_scale.x);
    const v = cam.vertical.scale((row + rand.float()) * cam.px_scale.y);
    const to_pixel = cam.top_left.add(u).subtract(v);

    const ray_dir: Vec3 = to_pixel.subtract(cam.position);
    return Ray{.origin = cam.position, .direction = ray_dir};
  }
};

//------------------------------------------------------------------------------
const HitInfo = struct
{
  point: Vec3,
  t: f32,
  surface_normal: Vec3,
  front_face: bool,
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
      .position = .{.x = 0.0, .y = -100.5, .z = -1.0},
      .radius = 100.0,
      .material =
      .{
        .Lambertian = .{.albedo = .{.x = 0.8, .y = 0.8, .z = 0.0}}
      }
    },
    .{
      .position = .{.x = 0.0, .y = 0.0, .z = -1.0},
      .radius = 0.5,
      .material =
      .{
        .Lambertian = .{.albedo = .{.x = 0.8, .y = 0.3, .z = 0.3}}
      }
    },
    .{
      .position = .{.x = 1.0, .y = 0.0, .z = -1.0},
      .radius = 0.5,
      .material =
      .{
        .Metal = .{.albedo = .{.x = 0.8, .y = 0.6, .z = 0.2}, .fuzz = 1.0}
      }
    },
    .{
      .position = .{.x = -1.0, .y = 0.0, .z = -1.0},
      .radius = 0.5,
      .material =
      .{
        .Dielectric = .{.refract_index = 1.5}
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

    const a: f32 = ray.direction.lengthSq();
    const half_b: f32 = Vec3.dot(displacement, ray.direction);
    const c: f32 = displacement.lengthSq() - (self.radius * self.radius);

    const discriminant: f32 = (half_b * half_b) - (a * c);
    if (discriminant < 0.0) { return null; }

    const discrim_root = math.sqrt(discriminant);
    const t1 = (-half_b - discrim_root) / a;
    const t2 = (-half_b + discrim_root) / a;

    if (closestValidT(t1, t2, t_min, t_max)) |t|
    {
      const point_at_t = ray.pointAtT(t);
      const outward_normal = point_at_t.subtract(self.position)
        .divideScalar(self.radius);

      return HitInfo
      {
        .point = point_at_t,
        .t = t,
        .surface_normal = outward_normal,
        .front_face = (outward_normal.dot(ray.direction) < 0.0),
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
  Dielectric,
};

const Material = union(MaterialTag)
{
  Lambertian: LambertianMaterial,
  Metal: MetalMaterial,
  Dielectric: DielectricMaterial,
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
    const normal_sphere_pos =  hit.point.add(hit.surface_normal);
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
  fuzz: f32,

  pub fn scatter(
    self: MetalMaterial,
    ray: Ray,
    hit: HitInfo,
    rand: *MyRand,
  ) ?ScatterInfo
  {
    const reflected_dir = ray.direction.reflect(hit.surface_normal);
    const fuzz_offset = rand.randomPointFromUnitSphere().scale(self.fuzz);
    const scattered_dir = reflected_dir.add(fuzz_offset);
    if (scattered_dir.dot(hit.surface_normal) > 0.0)
    {
      return ScatterInfo
      {
        .ray = Ray{.origin = hit.point, .direction = scattered_dir},
        .attenuation = self.albedo,
      };
    }

    return null;
  }
};

//------------------------------------------------------------------------------
const DielectricMaterial = struct
{
  refract_index: f32,

  pub fn schlick(cos_theta: f32, refract_index: f32) f32
  {
    const r = (1.0 - refract_index) / (1.0 + refract_index);
    const r0 = r * r;
    return r0 + (1.0 - r0) * math.pow(f32, (1.0 - cos_theta), 5.0);
  }

  pub fn scatter(
    self: DielectricMaterial,
    ray: Ray,
    hit: HitInfo,
    rand: *MyRand,
  ) ?ScatterInfo
  {
    // Adjust properties depending on whether ray is entering or exiting surface
    const ni_over_nt: f32 = if (hit.front_face)
      1.0 / self.refract_index
    else
      self.refract_index;

    const surface_normal = if (hit.front_face)
      hit.surface_normal
    else
      hit.surface_normal.scale(-1.0);

    // Obtain sin theta to check for total internal reflection
    const unit_v = ray.direction.normalized();
    const cos_theta = math.min(unit_v.scale(-1.0).dot(surface_normal), 1.0);
    const sin_theta = math.sqrt(1.0 - cos_theta*cos_theta);

    const will_reflect: bool = ( (ni_over_nt * sin_theta > 1.0)
      or (schlick(cos_theta, ni_over_nt) > rand.float()) );

    const scattered_dir: Vec3 = if (will_reflect)
      unit_v.reflect(surface_normal)
    else
      unit_v.refractTheta(cos_theta, surface_normal, ni_over_nt);

    return ScatterInfo
    {
      .ray = Ray{.origin = hit.point, .direction = scattered_dir},
      .attenuation = Vec3{.x = 1.0, .y = 1.0, .z = 1.0},
    };
  }
};

//------------------------------------------------------------------------------
pub fn calcColor(ray: Ray, scene: *const Scene, rand: *MyRand, call_depth: u32) Vec3
{
  if (call_depth > 50) { return Vec3{.x = 0.0, .y = 0.0, .z = 0.0}; }

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
      MaterialTag.Dielectric => |mat| mat.scatter(ray, hit, rand),
    };

    if (scatter_info) |scatter|
    {
      return calcColor(
        scatter.ray,
        scene,
        rand,
        call_depth + 1
      ).multiply(scatter.attenuation);
    }
    else
    {
      return Vec3{.x = 0.0, .y = 0.0, .z = 0.0};
    }
  }

  // Draw background
  const unit_direction = ray.direction.normalized();
  const t: f32 = (unit_direction.y + 1.0) * 0.5;

  const white_tone = Vec3.init(1.0, 1.0, 1.0).scale(1.0 - t);
  const blue_tone = Vec3.init(0.5, 0.7, 1.0).scale(t);
  return white_tone.add(blue_tone);
}

//------------------------------------------------------------------------------
pub fn main() anyerror!void
{
  // Scene to draw
  const scene_ptr = &basic_scene;

  // Output image details
  const image_filename = "test.bmp";
  const image_width: u32 =  1920;
  const image_height: u32 = 1080;

  // Rendering parameters
  const num_samples: u32 = 200;
  const num_samples_recip: f32 = 1.0 / @as(f32, num_samples);

  var rand = MyRand.create();

  const camera_pos = Vec3{.x = -2.0, .y = 2.0, .z = 1.0};
  const look_at_pos = Vec3{.x = 0.0, .y = 0.0, .z = -1.0};
  const vup = Vec3{.x = 0.0, .y = 1.0, .z = 0.0};
  const camera = Camera.create(
    camera_pos,
    look_at_pos,
    vup,
    image_width, image_height,
    20.0,
    @as(f32, image_width) / @as(f32, image_height),
  );

  // Output image
  var pixels: [image_width * image_height]Color = undefined;

  std.debug.warn("Reticulating splines...\n", .{});
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
