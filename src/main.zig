const c = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");

//------------------------------------------------------------------------------
const Color = packed struct {
    red: u8,
    green: u8,
    blue: u8,
    alpha: u8,

  pub const NUM_COMPONENTS: u32 = 4;
};

//------------------------------------------------------------------------------
pub fn main() anyerror!void {
  const IMAGE_FILENAME = "test.bmp";
  const IMAGE_WIDTH: u32 = 256;
  const IMAGE_HEIGHT: u32 = 256;

  var pixels: [IMAGE_WIDTH * IMAGE_HEIGHT]Color = undefined;

  // Fill with background color

  // Normalised range [0, 1]
  const x_range: f32 = @intToFloat(f32, IMAGE_WIDTH - 1);
  const y_range: f32 = @intToFloat(f32, IMAGE_HEIGHT - 1);

  // Scale pixel coordinate in range [0, 255]
  const x_scale: f32 = 255.0 / x_range;
  const y_scale: f32 = 255.0 / y_range;

  for (pixels) |*item, it|
  {
    const i = @intCast(u32, it);

    // The current pixel's coordinate position (row, col)
    const col: f32 = @intToFloat(f32, i % IMAGE_WIDTH);
    const row: f32 = @intToFloat(f32, i / IMAGE_WIDTH);

    item.red = @floatToInt(u8, col * x_scale);
    item.green = @floatToInt(u8, row * y_scale);
    item.blue  = 0x30;
    item.alpha = 0xFF;

    // std.debug.warn("Col:{}, Row:{}, X:{}, Y{}\n", .{col, row, x, y});
  }

  // Save the image to a file
  var f = c.stbi_write_bmp(
    IMAGE_FILENAME,
    @as(c_int, IMAGE_WIDTH),
    @as(c_int, IMAGE_HEIGHT),
    @as(c_int, Color.NUM_COMPONENTS),
    @ptrCast(*const c_void, &pixels[0])
  );

  std.debug.warn("All your codebase are belong to us.: {}\n", .{f});
}

//------------------------------------------------------------------------------
