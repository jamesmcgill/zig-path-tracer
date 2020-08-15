const c = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");

//------------------------------------------------------------------------------
pub fn main() anyerror!void {
  const IMAGE_FILENAME = "test.bmp";
  const IMAGE_WIDTH: u32 = 256;
  const IMAGE_HEIGHT: u32 = 256;
  const NUM_COMPONENTS: u32 = 4;

  var pixels: [IMAGE_WIDTH * IMAGE_HEIGHT * NUM_COMPONENTS]u8 = undefined;

  // Fill with background color
  // Byte-by-byte to avoid endianness concerns

  // Normalised range [0, 1]
  const xRange: f32 = @intToFloat(f32, IMAGE_WIDTH - 1);
  const yRange: f32 = @intToFloat(f32, IMAGE_HEIGHT - 1);

  // Scale pixel coordinate in range [0, 255]
  const xScale: f32 = 255.0 / xRange;
  const yScale: f32 = 255.0 / yRange;

  const componentsPerRow = NUM_COMPONENTS * IMAGE_WIDTH;
  for (pixels) |*item, it|
  {
    const i = @intCast(u32, it);
    // The current pixel's coordinate position (row, col)
    const row:f32 = @intToFloat(f32, i / componentsPerRow);
    const col:f32 = @intToFloat(f32, (i % componentsPerRow) / NUM_COMPONENTS);

    // std.debug.warn("Col:{}, Row:{}, X:{}, Y{}\n", .{col, row, x, y});
    const channel = i % NUM_COMPONENTS;
    item.* = switch (channel)
    {
      0 => @floatToInt(u8, col * xScale),  // Red
      1 => @floatToInt(u8, row * yScale),  // Green
      2 => 0x30,  // Blue
      3 => 0xff,  // Alpha
      else => 0x00,
    };
  }

  // Save the image to a file
  var f = c.stbi_write_bmp(
    IMAGE_FILENAME,
    @as(c_int, IMAGE_WIDTH),
    @as(c_int, IMAGE_HEIGHT),
    @as(c_int, NUM_COMPONENTS),
    @ptrCast(*const c_void, &pixels[0])
  );

  std.debug.warn("All your codebase are belong to us.: {}\n", .{f});
}

//------------------------------------------------------------------------------
