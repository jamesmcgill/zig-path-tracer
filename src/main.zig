const c = @cImport({
  @cInclude("stb_image_write.h");
});
const std = @import("std");

//------------------------------------------------------------------------------
pub fn main() anyerror!void {
  var raw: [256*256*4]u8 = undefined;
  for (raw) |*item, i|
  {
    const channel = @mod(i, 4);
    item.* = switch (channel)
    {
      0 => 0x00,  // R
      1 => 0x00,  // G
      2 => 0xff,  // B
      3 => 0xff,  // A
      else => 0x00,
    };
  }
  
  var f = c.stbi_write_bmp(
    "test.bmp",
    @as(c_int, 256),
    @as(c_int, 256),
    @as(c_int, 4),
    @ptrCast(*const c_void, &raw[0])
  );

  std.debug.warn("All your codebase are belong to us.: {}\n", .{f});
}

//------------------------------------------------------------------------------
