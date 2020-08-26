const c_import = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");
const Vec3 = @import("Vec3.zig").Vec3;
const Ray = @import("Ray.zig").Ray;

//------------------------------------------------------------------------------
const Color = packed struct {
    red: u8,
    green: u8,
    blue: u8,
    alpha: u8,

  pub const NUM_COMPONENTS: u32 = 4;
};

//------------------------------------------------------------------------------
pub fn hitsSphere(ray: Ray) bool {
  const sphere_pos = Vec3.init(0.0, 0.0, 110.0);
  const sphere_radius = 50.0;

  const displacement = ray.origin.subtract(sphere_pos);

  const a: f32 = Vec3.dot(ray.direction, ray.direction);
  const b: f32 = 2.0 * Vec3.dot(displacement, ray.direction);
  const c: f32 = Vec3.dot(displacement, displacement) -
    (sphere_radius * sphere_radius);

  const discriminant: f32 = (b * b) - (4.0 * a * c);
  return (discriminant >= 0.0);
}

//------------------------------------------------------------------------------
pub fn calcColor(ray: Ray) Vec3 {
  if (hitsSphere(ray))
  {
    return Vec3.init(1.0, 0.0, 0.0);
  }

  // Draw background
  const unit_direction = ray.direction.normalized();
  const t: f32 = (unit_direction.y + 1.0) * 0.5;

  const white_tone = Vec3.init(1.0, 1.0, 1.0).scale(t);
  const blue_tone = Vec3.init(0.5, 0.7, 1.0).scale(1.0 - t);
  return white_tone.add(blue_tone);
}

//------------------------------------------------------------------------------
pub fn main() anyerror!void {
  // Output image details
  const image_filename = "test.bmp";
  const image_width: u32 = 256;
  const image_height: u32 = 256;

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
    const i = @intCast(u32, it);

    // The current pixel's coordinate position (row, col)
    const col: f32 = @intToFloat(f32, i % image_width);
    const row: f32 = @intToFloat(f32, i / image_width);

    // The point that corresponds to on the frustum plane
    const to_x: f32 = x_offset + (col * x_scale);
    const to_y: f32 = y_offset + (row * y_scale);
    const to_pixel: Vec3 = Vec3.init(to_x, to_y, frustum_dist);

    const ray_dir: Vec3 = to_pixel.subtract(ray_origin);
    const ray = Ray.init(ray_origin, ray_dir);

    const color: Vec3 = calcColor(ray);
    item.red = @floatToInt(u8, color.x * 255.0);
    item.green = @floatToInt(u8, color.y * 255.0);
    item.blue = @floatToInt(u8, color.z * 255.0);
    item.alpha = 0xFF;

//    std.debug.warn("Col:{}, Row:{}, R:{}, G{}, B{}\n",
//      .{@floatToInt(u32, col), @floatToInt(u32, row),
//      item.red, item.green, item.blue});
  }

  // Save the image to a file
  var f = c_import.stbi_write_bmp(
    image_filename,
    @as(c_int, image_width),
    @as(c_int, image_height),
    @as(c_int, Color.NUM_COMPONENTS),
    @ptrCast(*const c_void, &pixels[0])
  );

  std.debug.warn("All your codebase are belong to us.: {}\n", .{f});
}

//------------------------------------------------------------------------------
