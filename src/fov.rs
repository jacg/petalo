/// The size and granularity of the Field of View (FOV) in which images should
/// be reconstructed

use crate::{Lengthf32, LengthU, LengthI, Pointf32};
use crate::{Length, Point, Vector, LOR, find_tof_peak, find_entry_point, voxel_size, first_boundaries};
use crate::index::{BoxDim_u, Index3_u, Index1_u, index1_to_3, index3_to_1};
use geometry::uom::mm_;
use geometry::{RatioPoint, in_base_unit};
use crate::guomc::ConstZero;

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
mod test_fov {
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
    fn test_voxel_centre(index: Index3_u, expected_position: [f32; 3]) {
        let fov = FOV::new((mm(4.0), mm(4.0), mm(4.0)), (2,2,2));
        let c = fov.voxel_centre(index);
        let c = [mm_(c.x), mm_(c.y), mm_(c.z)];
        assert_float_eq!(c, expected_position, ulps <= [1, 1, 1]);
    }
}

// --- Truncate float-based Lengthf32 to usize-based Lengthf32 --------------------------
#[inline(always)]
fn uom_floor(value: Length) -> LengthU { in_base_unit!(value.value.floor() as usize) }

#[inline(always)]
fn floor_f32(x: f32) -> usize { x.floor() as usize }

// --- Convert usize-based Lengthf32 to i32-based Lengthf32 -----------------------------
#[inline(always)]
fn uom_signed(value: LengthU) -> LengthI { in_base_unit!(value.value as i32) }

#[inline(always)]
fn signed_i32(x: usize) -> i32 { x as i32 }

/// Calculate information needed to keep track of progress across FOV:
/// voxel index and distance remaining until leaving the box
#[inline]
#[allow(clippy::identity_op)]
fn index_trackers(entry_point: RatioPoint, flipped: [bool; 3], [nx, ny, nz]: BoxDim_u) -> IndexTrackers {
    let entry_point: Pointf32 = entry_point.into();
    //use geometry::uom::uomcrate::ConstOne;
    //let one = ONE;
    let one = 1;

    // Find N-dimensional index of voxel at entry point.
    let [ix, iy, iz] = [floor_f32(entry_point.x),
                        floor_f32(entry_point.y),
                        floor_f32(entry_point.z)];

    // index is unsigned, but need signed values for delta_index
    let [ix, iy, iz] = [signed_i32(ix), signed_i32(iy), signed_i32(iz)];
    let [nx, ny, nz] = [signed_i32(nx), signed_i32(ny), signed_i32(nz)];

    // How much the 1d index changes along each dimension
    let delta_index = [
        1       * if flipped[0] { -one } else { one },
        nx      * if flipped[1] { -one } else { one },
        nx * ny * if flipped[2] { -one } else { one },
    ];

    // How many voxels remain before leaving FOV in each dimension
    let remaining = [
        nx - ix ,
        ny - iy ,
        nz - iz ,
    ];

    // 1d index into the 3d arrangement of voxels
    let [ix, iy, iz] = [
        if flipped[0] { nx - one - ix } else { ix },
        if flipped[1] { ny - one - iy } else { iy },
        if flipped[2] { nz - one - iz } else { iz },
    ];
    let index = index3_to_1([ix, iy, iz], [nx, ny, nz]);

    IndexTrackers { index, delta_index, remaining }
}

struct IndexTrackers {
    /// Current 1D index into 3D array of voxels
    index: i32,

    /// How the index changes along each dimension
    delta_index: [i32; 3],

    /// Voxels until edge of FOV in each dimension
    remaining: [i32; 3],
}

/// Flip axes to ensure that direction from p1 to p2 is non-negative in all
/// dimensions. Return p1 & p2 in flipped coordinate system, along with
/// knowledge of which axes were flipped, so that the indices of subsequently
/// found active voxels can be flipped back into the original coordinate system.
#[inline]
fn flip_axes(mut p1: Point, mut p2: Point) -> (Point, Point, [bool; 3]) {
    let zero = Length::ZERO;

    let dimensions = 3;
    let original_lor_direction: Vector = p2 - p1;
    let mut flipped = [false; 3];
    let mut flip_if_necessary = |n| {
        if original_lor_direction[n] < zero {
            p1[n] = - p1[n];
            p2[n] = - p2[n];
            flipped[n] = true;
        }
    };
    for d in 0..dimensions {
        flip_if_necessary(d);
    }
    (p1, p2, flipped)
}

/// Information about where the LOR enters the FOV and how to track the LOR's
/// traversal of the FOV.
pub struct FovHit {

    /// How far is the next voxel boundary in each direction.
    pub next_boundary: Vector,

    /// Voxel size expressed in LOR distance units: how far we must move along
    /// LOR to cross one voxel in any given dimension. Will be infinite for any
    /// axis which is parallel to the LOR.
    pub voxel_size   : Vector,

    /// 1D index of first voxel entered by the LOR.
    pub index        :  i32,

    /// Difference in 1D index between adjacent voxels, in each direction.
    pub delta_index  : [i32; 3],

    /// Number of voxels to be traversed along LOR, in each direction, before
    /// exiting FOV.
    pub remaining    : [i32; 3],

    /// Distance to the peak of the TOF gaussian.
    pub tof_peak     : Length,
}

/// Figure out if the LOR hits the FOV at all. If it does, calculate values
/// needed by `system_matrix_elements`.
#[inline]
pub fn lor_fov_hit(lor: &LOR, fov: FOV) -> Option<FovHit> {

    // Simplify expression of the algorithm by flipping axes so that the
    // direction from p1 to p2 is non-negative along all axes. Remember
    // which directions have been flipped, to recover correct voxel indices.
    let (p1, p2, flipped) = flip_axes(lor.p1, lor.p2);

    // If and where LOR enters FOV.
    let entry_point: Point = match fov.entry(p1, p2) {
        // If LOR misses the box, immediately return
        None => return None,
        // Otherwise, unwrap the point and continue
        Some(point) => point,
    };

    // How far the entry point is from the TOF peak
    let tof_peak = find_tof_peak(entry_point, p1, p2, lor.dt);

    // Express entry point in voxel coordinates: floor(position) = index of voxel.
    let entry_point: RatioPoint = find_entry_point(entry_point, fov);

    // Bookkeeping information needed during traversal of FOV
    let IndexTrackers {
        index,       // current 1d index into 3d array of voxels
        delta_index, // how the index changes along each dimension
        remaining,   // voxels until edge of FOV in each dimension
    } = index_trackers(entry_point, flipped, fov.n);

    // Voxel size expressed in LOR distance units: how far we must move along
    // LOR to cross one voxel in any given dimension. Will be infinite for any
    // axis which is parallel to the LOR.
    let voxel_size = voxel_size(fov, p1, p2);

    // At what position along LOR is the next voxel boundary, in any dimension.
    let next_boundary = first_boundaries(entry_point, voxel_size);

    // Return the values needed by `system_matrix_elements`
    let tof_peak = tof_peak;
    Some(FovHit { next_boundary, voxel_size, index, delta_index, remaining, tof_peak } )
}
