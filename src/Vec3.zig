const math = @import("std").math;

pub const Vec3 = struct {
    x: f32,
    y: f32,
    z: f32,

    //--------------------------------------------------------------------------
    pub fn init(x: f32, y: f32, z: f32) Vec3 {
        return Vec3{
            .x = x,
            .y = y,
            .z = z,
        };
    }

    //--------------------------------------------------------------------------
    pub fn add(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x + other.x,
            .y = self.y + other.y,
            .z = self.z + other.z,
        };
    }

    //--------------------------------------------------------------------------
    pub fn subtract(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x - other.x,
            .y = self.y - other.y,
            .z = self.z - other.z,
        };
    }

    //--------------------------------------------------------------------------
    pub fn scale(self: Vec3, scalar: f32) Vec3 {
        return Vec3{
            .x = self.x * scalar,
            .y = self.y * scalar,
            .z = self.z * scalar,
        };
    }

    //--------------------------------------------------------------------------
    pub fn lengthSq(self: Vec3) f32 {
        return self.dot(self);
    }

    //--------------------------------------------------------------------------
    pub fn length(self: Vec3) f32 {
        return math.sqrt(self.lengthSq());
    }

    //--------------------------------------------------------------------------
    pub fn normalized(self: Vec3) Vec3 {
        return self.scale(1.0 / self.length());
    }

    //--------------------------------------------------------------------------
    pub fn dot(self: Vec3, other: Vec3) f32 {
        return self.x * other.x + self.y * other.y + self.z * other.z;
    }

    //--------------------------------------------------------------------------
};
