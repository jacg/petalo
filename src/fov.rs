/// The size and granularity of the Field of View (FOV) in which images should
/// be reconstructed

use crate::types::Lengthf32;
use crate::types::{Length, Point, Vector};
use crate::index::{BoxDim_u, Index3_u, Index1_u, index1_to_3};
use geometry::uom::mm_;

#[derive(Clone, Copy, Debug)]
pub struct FOV {
    pub half_width: Vector,
    pub n: BoxDim_u,
    pub voxel_size: Vector,
}

impl FOV {

    pub fn new(
        full_size: (Length, Length, Length),
        (nx, ny, nz): (usize, usize, usize)
    ) -> Self {
        let (dx, dy, dz) = full_size;
        let half_width = Vector::new(dx/2.0, dy/2.0, dz/2.0);
        let n = [nx, ny, nz];
        let voxel_size = Self::voxel_size(n, half_width);
            Self { half_width, n, voxel_size }
    }

    fn voxel_size(n: BoxDim_u, half_width: Vector) -> Vector {
        let full_width = half_width * 2.0;
        Vector::new(full_width[0] / n[0] as f32,
                    full_width[1] / n[1] as f32,
                    full_width[2] / n[2] as f32,
        )
    }

    /// Find centre of voxel with given 3D index
    pub fn voxel_centre(&self, i: Index3_u) -> Point {
        //i.map(|n| n as f64 + 0.5).component_mul(&self.voxel_size).into()
        let s = self.voxel_size;
        Point::new((i[0] as Lengthf32 + 0.5) * s.x - self.half_width[0],
                   (i[1] as Lengthf32 + 0.5) * s.y - self.half_width[1],
                   (i[2] as Lengthf32 + 0.5) * s.z - self.half_width[2],)
    }

    /// Find centre of voxel with given 1D index
    pub fn voxel_centre1(&self, i: Index1_u) -> Point {
        self.voxel_centre(index1_to_3(i, self.n))
    }

    pub fn entry(&self, p1: Point, p2: Point) -> Option<Point> {

        use ncollide3d::query::RayCast;
        use ncollide3d::shape::Cuboid;

        type Ray      = ncollide3d::query::Ray    <Lengthf32>;
        type Isometry = ncollide3d::math::Isometry<Lengthf32>;

        let lor_direction = (p2 - p1).normalize();
        let lor_length    = (p2 - p1).norm();
        let lor: Ray = Ray::new(p1.into(), lor_direction.into());
        let iso: Isometry = Isometry::identity();
        Cuboid::new(self.half_width.into())
            .toi_with_ray(&iso, &lor, mm_(lor_length), true)
            .map(|toi| lor.origin + lor.dir * toi)
            .map(Into::into)
    }

}

#[cfg(test)]
mod test_voxel_box {
    use super::*;
    use rstest::rstest;
    use geometry::uom::mm;
    use float_eq::assert_float_eq;

    #[rstest(/**/ index,   expected_position,
             case([0,0,0], [-1.0, -1.0, -1.0]),
             case([0,0,1], [-1.0, -1.0,  1.0]),
             case([0,1,0], [-1.0,  1.0, -1.0]),
             case([0,1,1], [-1.0,  1.0,  1.0]),
             case([1,0,0], [ 1.0, -1.0, -1.0]),
             case([1,0,1], [ 1.0, -1.0,  1.0]),
             case([1,1,0], [ 1.0,  1.0, -1.0]),
             case([1,1,1], [ 1.0,  1.0,  1.0]),
    )]
    fn test_voxel_centre(index: Index3_u, expected_position: [Lengthf32; 3]) {
        let fov = FOV::new((mm(4.0), mm(4.0), mm(4.0)), (2,2,2));
        let c = fov.voxel_centre(index);
        let c = [mm_(c.x), mm_(c.y), mm_(c.z)];
        assert_float_eq!(c, expected_position, ulps <= [1, 1, 1]);
    }
}
