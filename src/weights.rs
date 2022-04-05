//! Find the weights and indices of the active voxels along a single Line Of
//! Response LOR.
//!
//! The algorithm is centred around two key simplifications:
//!
//! 1. Express the voxel size in terms of the components of the LOR's direction
//!    vector. This allows trivial calculation of how far we must move along the
//!    LOR before reaching a voxel boundary, in any dimension.
//!
//! 2. Exploit symmetry to simplify dealing with directions: flip axes so that
//!    the direction of the LOR has non-negative components. The algorithm can
//!    then assume that all progress is in the positive direction. Any voxels
//!    indices calculated by the algorithm, must be flipped back to the original
//!    coordinate system.

use geometry::in_base_unit;
use crate::types::{BoxDim_u, Index1_u, Index3_u, Index3Weightf32, Lengthf32, Pointf32, Vectorf32, ns_to_mm};
use crate::types::{Length, LengthU, LengthI, PerLength, Time, C,
                   Point, Vector, Ratio, RatioPoint, RatioVec};

use geometry::uom::{mm, mm_, ns_};
use crate::gauss::make_gauss_option;
use crate::mlem::{index3_to_1, index1_to_3};

// ------------------------------ TESTS ------------------------------
#[cfg(test)]
use float_eq::assert_float_eq;

#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use assert_approx_eq::assert_approx_eq;
    use crate::types::TWOPIf32;
    use geometry::uom::ratio;

    // --------------------------------------------------------------------------------
    // This set of hand-picked values should be easy to verify by humans. The
    // test performs two checks:
    //
    // 1. The sum of the LOR-lengths within individual voxels equals the
    //    expected total length of LOR in the whole FOV.
    //
    // 2. The indices of the voxels traversed by the LOR are as expected.
    #[rstest(/**/      p1       ,      p2      ,    size     ,  n   ,  length  , expected_voxels,
             // symmetric 3x3, diagonal LOR under all four axis flip combinations
             case((-30.0, -30.0), ( 30.0, 30.0), (10.0, 10.0), (3,3), 14.142135, vec![(0,0), (1,1), (2,2)]),
             case(( 30.0, -30.0), (-30.0, 30.0), (10.0, 10.0), (3,3), 14.142135, vec![(2,0), (1,1), (0,2)]),
             case((-30.0,  30.0), ( 30.0,-30.0), (10.0, 10.0), (3,3), 14.142135, vec![(0,2), (1,1), (2,0)]),
             case(( 30.0,  30.0), (-30.0,-30.0), (10.0, 10.0), (3,3), 14.142135, vec![(2,2), (1,1), (0,0)]),
             // like case 1, but with asymmetric voxels
             case((-30.0, -30.0), ( 30.0, 30.0), (10.0, 10.0), (3,2), 14.142135, vec![(0,0), (1,0), (1,1), (2,1)]),
             case((-30.0, -30.0), ( 30.0, 30.0), (10.0, 10.0), (2,3), 14.142135, vec![(0,0), (0,1), (1,1), (1,2)]),
             // vertical / horizontal off-centre LOR
             case((  5.4, -20.0), (  5.4, 10.0), (11.0,  9.0), (9,4),  9.0     , vec![(8,0), (8,1), (8,2), (8,3)]),
             case((-15.0,  -4.0), ( 15.0, -4.0), ( 8.0, 10.0), (4,3),  8.0     , vec![(0,0), (1,0), (2,0), (3,0)]),
    )]
    fn hand_picked(p1:   (Lengthf32, Lengthf32),
                   p2:   (Lengthf32, Lengthf32),
                   size: (Lengthf32, Lengthf32),
                   n: (usize, usize),
                   length: Lengthf32,
                   expected_voxels: Vec<(usize, usize)>) {

        let p1 = (mm(p1.0), mm(p1.1));
        let p2 = (mm(p2.0), mm(p2.1));

        let p1 = Point::new(p1.0, p1.1, mm(0.0));
        let p2 = Point::new(p2.0, p2.1, mm(0.0));
        let fov = FOV::new((mm(size.0), mm(size.1), mm(1.0)), (n.0, n.1, 1));

        // Values to plug in to visualizer:
        let lor = LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0));
        let command = crate::visualize::vislor_command(&fov, &lor);
        println!("\nTo visualize this case, run:\n{}\n", command);

        // Collect hits
        let hits: Vec<Index3Weightf32> = LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)).active_voxels(&fov, None, None);

        // Diagnostic output
        for (is, l) in &hits { println!("  ({} {})   {}", is[0], is[1], l) }

        // Check total length through FOV
        let total_length: Lengthf32 = hits.iter()
            .map(|(_index, weight)| weight)
            .sum();
        assert_approx_eq!(total_length, length);

        // Check voxels hit
        let voxels: Vec<(usize, usize)> = hits.into_iter()
            .map(|(index, _weight)| (index[0], index[1]))
            .collect();
        assert_eq!(voxels, expected_voxels)
    }

    // --------------------------------------------------------------------------------
    use proptest::prelude::*;
    // This property-based test generates random test cases and verifies that
    // the total length of the LOR in the FOV equals the sum of its lengths in
    // the individual voxels.
    proptest! {
        #[test]
        fn sum_of_weights_equals_length_through_box(
            // Activated sensor positions
            r        in  200.0..(300.0 as Lengthf32),
            p1_angle in 0.0..(1.0 as Lengthf32), // around the circle
            p2_delta in 0.1..(0.9 as Lengthf32), // relative to p1_angle
            p1_z     in -200.0..(200.0 as Lengthf32),
            p2_z     in -200.0..(200.0 as Lengthf32),
            // Field of View
            dx in  100.0..(150.0 as Lengthf32),
            dy in  100.0..(150.0 as Lengthf32),
            dz in  100.0..(190.0 as Lengthf32),
            nx in  5..50_usize,
            ny in  5..50_usize,
            nz in  5..90_usize,
        ) {
            let p1_theta: Lengthf32 = p1_angle * TWOPIf32;
            let p2_theta: Lengthf32 = p1_theta + (p2_delta * TWOPIf32);
            let p1 = Pointf32::new(r * p1_theta.cos(), r * p1_theta.sin(), p1_z);
            let p2 = Pointf32::new(r * p2_theta.cos(), r * p2_theta.sin(), p2_z);
            let fov = FOV::new((mm(dx), mm(dy), mm(dz)), (nx, ny, nz));

            // Values to plug in to visualizer:
            let lor = LOR::new(Time::ZERO, Time::ZERO, p1.into(), p2.into(), ratio(1.0));
            let command = crate::visualize::vislor_command(&fov, &lor);
            println!("\nTo visualize this case, run:\n{}\n", command);

            let summed: Lengthf32 = LOR::new(Time::ZERO, Time::ZERO, p1.into(), p2.into(), ratio(1.0))
                .active_voxels(&fov, None, None)
                .into_iter()
                .inspect(|(i, l)| println!("  ({} {} {}) {}", i[0], i[1], i[2], l))
                .map(|(_index, weight)| weight)
                .sum();

            let a = fov.entry(&p1.into(), &p2.into());
            let b = fov.entry(&p2.into(), &p1.into());

            let in_one_go = match (a,b) {
                (Some(a), Some(b)) => (a - b).magnitude(),
                _ => mm(0.0)
            };

            assert_float_eq!(summed, mm_(in_one_go), rel <= 1e-3);

        }
    }
}

// ---------------------- Implementation -----------------------------------------

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
    let entry_point: Point = match fov.entry(&p1, &p2) {
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

/// For a single LOR, place the weights and indices of the active voxels in the
/// `weights` and `indices` parameters. Using output parameters rather than
/// return values, because this function is called in the inner loop, and
/// allocating the vectors of results repeatedly, had a noticeable impact on
/// performance.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn system_matrix_elements(
    indices: &mut Vec<usize>,
    weights: &mut Vec<Lengthf32>,
    mut next_boundary: Vectorf32,
    voxel_size: Vector,
    mut index: i32,
    delta_index: [i32; 3],
    mut remaining: [i32; 3],
    tof_peak: Length,
    tof: &Option<impl Fn(Length) -> PerLength>) {

    // How far we have moved since entering the FOV
    let mut here = F32_LENGTH_ZERO;

    loop {
        // Which voxel boundary will be hit next, and its position
        let (dimension, boundary_position) = next_boundary.argmin();

        // The weight is the length of LOR in this voxel
        let mut weight = boundary_position - here;

        // If TOF enabled, adjust weight
        if let Some(gauss) = &tof {
            let g: PerLength = gauss(mm(here) - tof_peak);
            // Turn into dimensionless number: TODO normalization
            let g: f32 = (mm(1000.0) * g).get::<geometry::uom::uomcrate::si::ratio::ratio>();
            weight *= g;
        }

        // Store the index and weight of the voxel we have just crossed
        if weight > 0.0 {
            indices.push(index as usize);
            weights.push(weight);
        }

        // Move along LOR until it leaves this voxel
        here = boundary_position;

        // Find the next boundary in this dimension
        next_boundary[dimension] += Vectorf32::from(voxel_size)[dimension];

        // Move index across the boundary we are crossing
        index += delta_index[dimension];
        remaining[dimension] -= 1;

        // If we have traversed the whole FOV, we're finished
        if remaining[dimension] == 0 { break; }
    }
}

use crate::types::guomc::ConstZero;
const F32_LENGTH_ZERO: Lengthf32 = 0.0;

const EPS: Ratio = in_base_unit!(1e-5);

/// The point at which the LOR enters the FOV, expressed in a coordinate
/// system with one corner of the FOV at the origin.
#[inline]
fn find_entry_point(entry_point: Point, fov: FOV) -> RatioPoint {
    // Transform coordinates to align box with axes: making the lower boundaries
    // of the box lie on the zero-planes.
    (entry_point + fov.half_width)

        // Express entry point in voxel coordinates: floor(position) = index of
        // voxel.
        .component_div(fov.voxel_size)

        // Floating-point subtractions which should give zero, usually miss very
        // slightly: if this error is negative, the next step (which uses floor)
        // will pick the wrong voxel. Work around this problem by assuming that
        // anything very close to zero is exactly zero.
        .map(|x| if x.abs() < EPS { Ratio::ZERO } else { x })

    // Express entry point in voxel coordinates: floor(position) = index of voxel.
    // TODO: figure out if we should support Point * Vector -> Point  (affine * vector -> affine)
    // NOTE: this should be Point<Ratio> rater than Point<Lengthf32>

    // let mut entry_point = Pointf32::new(
    //     entry_point[0] / voxel_size[0],
    //     entry_point[1] / voxel_size[1],
    //     entry_point[2] / voxel_size[2],
    // );

    // // Floating-point subtractions which should give zero, usually miss very
    // // slightly: if this error is negative, the next step (which uses floor)
    // // will pick the wrong voxel. Work around this problem by assuming that
    // // anything very close to zero is exactly zero.
    // entry_point.iter_mut().for_each(|x| if x.abs() < EPSf32 { *x = 0.0 });
    // entry_point.into()

}


/// Distance from entry point to the LOR's TOF peak
#[inline]
fn find_tof_peak(entry_point: Point, p1: Point, p2: Point, dt: Time) -> Length {
    let half_lor_length = (p1 - p2).norm() / 2.0;
    let tof_shift = C * dt / 2.0; // NOTE ignoring refractive index
    let p1_to_peak = half_lor_length - tof_shift;
    let p1_to_entry = (entry_point - p1).norm();
    p1_to_peak - p1_to_entry
}


/// Distances from entry point to the next voxel boundaries, in each dimension
#[inline]
fn first_boundaries(entry_point: RatioPoint, voxel_size: Vector) -> Vector {
    use geometry::uom::uomcrate::si::ratio::ratio;
    // How far have we penetrated into this voxel, along any axis
    let frac_done: RatioVec = entry_point - entry_point.map(|x| x.floor::<ratio>());
    // Distances remaining to the nearest boundaries
    (RatioVec::new(1.0, 1.0, 1.0) - frac_done) * voxel_size
}

/// Voxel size expressed in LOR distance units: how far we must move along LOR
/// to cross one voxel in any given dimension. Will be infinite for any axis
/// which is parallel to the LOR.
#[inline]
fn voxel_size(fov: FOV, p1: Point, p2: Point) -> Vector {
    // TODO: The units are a bit dodgy here. See the TODOs for
    // Vector::{normalize,component_div}
    let lor_direction = (p2-p1).normalize();
    fov.voxel_size.component_div(&lor_direction)
}

use geometry::Quantity;

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


//--------------------------------------------------------------------------------
/// The size and granularity of the region in which images should be
/// reconstructed
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

    pub fn entry(&self, p1: &Point, p2: &Point) -> Option<Point> {

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

//--------------------------------------------------------------------------------
/// Line Of Response.
///
/// 2 spacetime vectors indicating the positions and times of coincident
/// detector element activations
#[derive(Clone, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct LOR {
    pub p1: Point,
    pub p2: Point,
    pub dt: Time,
    /// Scatter and random corrections, which appear as an additive contribution
    /// to the sinogram bin in the MLEM forward projection. In order to make it
    /// compatible with a single LOR (rather than many LORs in a sinogram bin)
    /// it is expressed here as a *multiplicative* factor.
    pub additive_correction: Ratio,
}

impl LOR {
    pub fn new(t1: Time, t2: Time, p1: Point, p2: Point, additive_correction: Ratio) -> Self {
        Self { p1, p2, dt: t2 - t1, additive_correction }
    }

    pub fn from_components((t1, t2): (Time, Time),
                           (x1, y1, z1): (Length, Length, Length),
                           (x2, y2, z2): (Length, Length, Length),
                           additive_correction: Ratio
                          ) -> Self
    {
        Self::new(t1, t2, Point::new(x1,y1,z1), Point::new(x2,y2,z2), additive_correction)
    }

    pub fn active_voxels(&self, fov: &FOV, cutoff: Option<Ratio>, sigma: Option<Time>) -> Vec<Index3Weightf32> {
        let tof = make_gauss_option(sigma, cutoff);
        let mut weights = vec![];
        let mut indices = vec![];
        match lor_fov_hit(self, *fov) {
            None => (),
            Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {
                system_matrix_elements(
                    &mut indices, &mut weights,
                    next_boundary.into(), voxel_size,
                    index, delta_index, remaining,
                    tof_peak, &tof
                );

            }
        }
        indices.into_iter()
            .map(|i| index1_to_3(i, fov.n))
            .zip(weights.into_iter())
            .collect()
    }
}

use core::fmt;
impl fmt::Display for LOR {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (p, q) = (self.p1, self.p2);
        write!(f, "<LOR ({:8.2} {:8.2} {:8.2}) ({:8.2} {:8.2} {:8.2}) {:7.2}ns {:7.2}mm /{:7.2} >",
               mm_(p.x), mm_(p.y), mm_(p.z),
               mm_(q.x), mm_(q.y), mm_(q.z),
               ns_(self.dt), ns_to_mm(ns_(self.dt)) / 2.0,
               mm_((p-q).norm())
        )
    }
}
