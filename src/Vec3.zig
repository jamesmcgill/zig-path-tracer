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
    pub fn multiply(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x * other.x,
            .y = self.y * other.y,
            .z = self.z * other.z,
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
    pub fn sqrt(self: Vec3) Vec3 {
        return Vec3{
            .x = math.sqrt(self.x),
            .y = math.sqrt(self.y),
            .z = math.sqrt(self.z),
        };
    }

    //--------------------------------------------------------------------------
    pub fn reflect(self: Vec3, surface_normal: Vec3) Vec3 {
        const proj = self.dot(surface_normal);
        return self.subtract(surface_normal.scale(proj).scale(2.0));
    }

    //--------------------------------------------------------------------------
    pub fn refract(self: Vec3, surface_normal: Vec3, ni_over_nt: f32) Vec3 {
        const unit_v = self.normalized().scale(-1.0);
        const cos_theta = unit_v.dot(surface_normal);

        const proj_n = surface_normal.scale(cos_theta);
        const r_out_perp = unit_v.add(proj_n).scale(ni_over_nt);

        const parallel_len = -math.sqrt(math.fabs(1.0 - r_out_perp.lengthSq()));
        const r_out_parallel = surface_normal.scale(parallel_len);

        return r_out_perp.add(r_out_parallel);
    }

    //--------------------------------------------------------------------------
};
