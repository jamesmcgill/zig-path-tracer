const Vec3 = @import("Vec3.zig").Vec3;

pub const Ray = struct {
    origin: Vec3,
    direction: Vec3,
    time: f32 = 0.0,

    //--------------------------------------------------------------------------
    pub fn init(origin: Vec3, direction: Vec3, time: f32) Ray {
        return Ray{
            .origin = origin,
            .direction = direction,
            .time = time,
        };
    }

    //--------------------------------------------------------------------------
    pub fn pointAtT(self: Ray, t: f32) Vec3 {
        return self.origin.add(self.direction.scale(t));
    }

    //--------------------------------------------------------------------------
};
