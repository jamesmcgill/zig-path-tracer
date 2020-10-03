const std = @import("std");
const math = @import("std").math;
const assert = @import("std").debug.assert;

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
    pub fn multiply(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x * other.x,
            .y = self.y * other.y,
            .z = self.z * other.z,
        };
    }

    //--------------------------------------------------------------------------
    pub fn divide(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x / other.x,
            .y = self.y / other.y,
            .z = self.z / other.z,
        };
    }

    //--------------------------------------------------------------------------
    pub fn addScalar(self: Vec3, val: f32) Vec3 {
        return Vec3{
            .x = self.x + val,
            .y = self.y + val,
            .z = self.z + val,
        };
    }

    //--------------------------------------------------------------------------
    pub fn subtractScalar(self: Vec3, val: f32) Vec3 {
        return Vec3{
            .x = self.x - val,
            .y = self.y - val,
            .z = self.z - val,
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
    pub fn divideScalar(self: Vec3, scalar: f32) Vec3 {
        return Vec3{
            .x = self.x / scalar,
            .y = self.y / scalar,
            .z = self.z / scalar,
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
    // Scale a normalized vector in range [-1, 1] to range [0, 1]
    //--------------------------------------------------------------------------
    // TODO: find a better name (requires a normal vector) eg positiveScaleNormal
    pub fn rescaled(self: Vec3) Vec3 {
        return self.addScalar(1.0).scale(0.5);
    }

    //--------------------------------------------------------------------------
    pub fn dot(self: Vec3, other: Vec3) f32 {
        return self.x * other.x + self.y * other.y + self.z * other.z;
    }

    //--------------------------------------------------------------------------
    pub fn cross(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = (self.y * other.z) - (self.z * other.y),
            .y = (self.z * other.x) - (self.x * other.z),
            .z = (self.x * other.y) - (self.y * other.x),
        };
    }

    //--------------------------------------------------------------------------
    pub fn sqrt(self: Vec3) Vec3 {
        return Vec3{
            .x = math.sqrt(self.x),
            .y = math.sqrt(self.y),
            .z = math.sqrt(self.z),
        };
    }

    //--------------------------------------------------------------------------
    pub fn reflect(self: Vec3, surface_normal: Vec3) Vec3 {
        const double_height = self.dot(surface_normal) * -2;
        return self.add(surface_normal.scale(double_height));
    }

    //--------------------------------------------------------------------------
    pub fn refract(self: Vec3, surface_normal: Vec3, ni_over_nt: f32) Vec3 {
        const unit_v = self.normalized();
        const cos_theta = math.min(unit_v.scale(-1.0).dot(surface_normal), 1.0);

        return unit_v.refractTheta(cos_theta, surface_normal, ni_over_nt);
    }

    //--------------------------------------------------------------------------
    pub fn refractTheta(self: Vec3, in_cos_theta: f32, surface_normal: Vec3, ni_over_nt: f32) Vec3 {
        assert(math.approxEq(f32, 1.0, self.lengthSq(), 0.000001)); // self must be normalized
        const in_ortho_len = surface_normal.scale(in_cos_theta);
        const ortho = self.add(in_ortho_len).scale(ni_over_nt);

        const parallel_len = -math.sqrt(math.fabs(1.0 - ortho.lengthSq()));
        const parallel = surface_normal.scale(parallel_len);

        return ortho.add(parallel);
    }

    //--------------------------------------------------------------------------
    pub fn print(self: Vec3) void {
        std.debug.warn("({d:.3},{d:.3},{d:.3})", .{ self.x, self.y, self.z });
    }

    //--------------------------------------------------------------------------
};
