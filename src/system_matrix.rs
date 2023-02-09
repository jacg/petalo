//! Find indices and weights of the voxels coupled to a single LOR
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

use units::{
    C, Length, Ratio, PerLength, Time,
    mm, mm_, ratio_,
    in_base_unit,
    todo::{Lengthf32, Weightf32},
};

use crate::{
    BoxDim_u, Index1_u, LOR, Point, Pointf32, RatioPoint, RatioVec, Vector,
    config::mlem::Tof,
    fov::FOV,
    gauss::{Gaussian, make_gauss_option},
    image::{Image, ImageData},
    index::index3_to_1,
};


// ------------------------------ TESTS ------------------------------
#[cfg(test)]
use float_eq::assert_float_eq;

#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use units::{TWOPI, ratio};
    use crate::index::index1_to_3;

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

        // Collect voxels traversed by LOR
        let hits = Siddon::new_system_matrix_row(&LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)), &fov, None);

        // Utility for converting 1D-index to 3D-index
        let as_3d = |i| index1_to_3(i, [n.0, n.1, 1]);

        // Diagnostic output
        for (i, l) in &hits { println!("  ({} {})   {l}", as_3d(i)[0], as_3d(i)[1]) }

        // Check total length through FOV
        let total_length: Lengthf32 = hits.iter()
            .map(|(_index, weight)| weight)
            .sum();
        assert_float_eq!(total_length, length, ulps <= 1);

        // Check voxels hit
        let voxels: Vec<(usize, usize)> = hits.into_iter()
            .map(|(j, _weight)| (as_3d(j)[0], as_3d(j)[1]))
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
            let (r, p1_z, p2_z) = (mm(r), mm(p1_z), mm(p2_z));
            let p1_theta = p1_angle * TWOPI;
            let p2_theta = p1_theta + (p2_delta * TWOPI);
            let p1 = Point::new(r * p1_theta.cos(), r * p1_theta.sin(), p1_z);
            let p2 = Point::new(r * p2_theta.cos(), r * p2_theta.sin(), p2_z);
            let fov = FOV::new((mm(dx), mm(dy), mm(dz)), (nx, ny, nz));

            // Convert 1D-index to 3D-index
            let as_3d = |i| index1_to_3(i, [nx, ny, nz]);

            // Values to plug in to visualizer:
            let lor = LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0));
            let command = crate::visualize::vislor_command(&fov, &lor);
            println!("\nTo visualize this case, run:\n{}\n", command);

            let summed: Lengthf32 = Siddon::new_system_matrix_row(
                &LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)), &fov, None
            ).into_iter()
             .inspect(|&(i, l)| println!("  ({} {} {}) {}", as_3d(i)[0], as_3d(i)[1], as_3d(i)[2], l))
             .map(|(_index, weight)| weight)
             .sum();

            let a = fov.entry(p1, p2);
            let b = fov.entry(p2, p1);

            let in_one_go = match (a,b) {
                (Some(a), Some(b)) => (a - b).magnitude(),
                _ => mm(0.0)
            };

            assert_float_eq!(summed, mm_(in_one_go), rel <= 1e-3);

        }
    }
}

// ---------------------- Implementation -----------------------------------------
pub type SystemMatrixElement = (Index1_u, Weightf32);

pub struct SystemMatrixRow(pub Vec<SystemMatrixElement>);

impl SystemMatrixRow {
    pub fn iter(&self) -> std::slice::Iter<SystemMatrixElement> { self.0.iter() }
    pub fn clear(&mut self) { self.0.clear(); }
}

impl IntoIterator for SystemMatrixRow {
    type Item = SystemMatrixElement;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a SystemMatrixRow {
    type Item = SystemMatrixElement;
    type IntoIter = std::iter::Cloned<std::slice::Iter<'a, Self::Item>>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().cloned()
    }
}

// ---- Projector Trait ----------------------------------------------------------------------
pub trait Projector {
    fn project_one_lor<'i, 'g>(fold_state: FoldState<'i, 'g>, lor: &LOR) -> FoldState<'i, 'g>;

    // Sparse storage of the slice through the system matrix which corresponds
    // to the current LOR. Allocating these anew for each LOR had a noticeable
    // runtime cost, so we create them up-front and reuse them.
    // This should probably have a default implementation
    fn buffers(fov: FOV) -> SystemMatrixRow;
}

// ---- Siddon Projector ----------------------------------------------------------------------
pub struct Siddon;

impl Projector for Siddon {

    fn project_one_lor<'i, 'g>(fold_state: FoldState<'i, 'g>, lor: &LOR) -> FoldState<'i, 'g> {
        project_one_lor(fold_state, lor)
    }

    fn buffers(fov: FOV) -> SystemMatrixRow {
        let [nx, ny, nz] = fov.n;
        let max_number_of_coupled_voxels_possible = nx + ny + nz - 2;
        SystemMatrixRow(Vec::with_capacity(max_number_of_coupled_voxels_possible))
    }

}

impl Siddon {
    // TODO  FOV and TOF should become construction-time arguments
    pub fn new_system_matrix_row(lor: &LOR, fov: &FOV, tof: Option<Tof>) -> SystemMatrixRow {
        let tof = make_gauss_option(tof);
        let mut system_matrix_row = Self::buffers(*fov);
        match lor_fov_hit(lor, *fov) {
            None => (),
            Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {
                Self::update_smatrix_row(
                    &mut system_matrix_row,
                    next_boundary, voxel_size,
                    index, delta_index, remaining,
                    tof_peak, &tof
                );
            }
        }
        system_matrix_row
    }

    /// For a single LOR, place the weights and indices of the coupled voxels in
    /// `system_matrix_row` parameter. Using output parameters rather than return
    /// values, because this function is called in the inner loop, and allocating
    /// the vectors of results repeatedly, had a noticeable impact on performance.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    pub fn update_smatrix_row(
        system_matrix_row: &mut SystemMatrixRow,
        mut next_boundary: Vector,
        voxel_size: Vector,
        mut index: i32,
        delta_index: [i32; 3],
        mut remaining: [i32; 3],
        tof_peak: Length,
        tof: &Option<Gaussian>
    ) {
        // Throw away previous LOR's values
        system_matrix_row.clear();

        // How far we have moved since entering the FOV
        let mut here = Length::ZERO;

        loop {
            // Which voxel boundary will be hit next, and its position
            let (dimension, boundary_position) = next_boundary.argmin();

            // The weight is the length of LOR in this voxel
            let mut weight = boundary_position - here;

            // If TOF enabled, adjust weight
            if let Some(gauss) = &tof {
                let g: PerLength = gauss.call(here - tof_peak);
                // TODO Normalization
                let completely_arbitrary_factor = 666.0;
                let g: f32 = ratio_(mm(completely_arbitrary_factor) * g);
                weight *= g;
            }

            // Store the index and weight of the voxel we have just crossed
            if weight > Length::ZERO {
                system_matrix_row.0.push(((index as usize), (mm_(weight))));
            }

            // Move along LOR until it leaves this voxel
            here = boundary_position;

            // Find the next boundary in this dimension
            next_boundary[dimension] += voxel_size[dimension];

            // Move index across the boundary we are crossing
            index += delta_index[dimension];
            remaining[dimension] -= 1;

            // If we have traversed the whole FOV, we're finished
            if remaining[dimension] == 0 { break; }
        }
    }

}


use units::uom::ConstZero;

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
    use units::uom::si::ratio::ratio;
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
    fov.voxel_size.component_div(lor_direction)
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

/// Calculate information needed to keep track of progress across FOV:
/// voxel index and distance remaining until leaving the box
#[inline]
#[allow(clippy::identity_op)]
fn index_trackers(entry_point: RatioPoint, flipped: [bool; 3], [nx, ny, nz]: BoxDim_u) -> IndexTrackers {
    let entry_point: Pointf32 = entry_point.into();
    //use units::uom::ConstOne;
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

#[inline(always)]
fn floor_f32(x: f32) -> usize { x.floor() as usize }

#[inline(always)]
fn signed_i32(x: usize) -> i32 { x as i32 }

/// Flip axes to ensure that direction from p1 to p2 is non-negative in all
/// dimensions. Return p1 & p2 in flipped coordinate system, along with
/// knowledge of which axes were flipped, so that the indices of subsequently
/// found voxels can be flipped back into the original coordinate system.
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

pub type FoldState<'i, 'g> = (ImageData, SystemMatrixRow , &'i Image, &'g Option<Gaussian>);

pub fn project_one_lor<'i, 'g>(state: FoldState<'i, 'g>, lor: &LOR) -> FoldState<'i, 'g>{
    let (mut backprojection, mut system_matrix_row, image, tof) = state;

    // Need to return the state from various match arms
    macro_rules! return_state { () => (return (backprojection, system_matrix_row, image, tof)); }

    // Analyse point where LOR hits FOV
    match lor_fov_hit(lor, image.fov) {

        // LOR missed FOV: nothing to be done
        None => return_state!(),

        // Data needed by `system_matrix_elements`
        Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

            // Find non-zero elements (voxels coupled to this LOR) and their values
            Siddon::update_smatrix_row(
                &mut system_matrix_row,
                next_boundary, voxel_size,
                index, delta_index, remaining,
                tof_peak, &tof
            );

            // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
            for (i, _) in &system_matrix_row {
                if i >= backprojection.len() { return_state!(); }
            }

            // Forward projection of current image into this LOR
            let projection = forward_project(&system_matrix_row, image) * lor.additive_correction;

            // Backprojection of LOR onto image
            back_project(&mut backprojection, &system_matrix_row, ratio_(projection));
            return_state!();
        }
    }
}

#[inline]
pub fn forward_project(system_matrix_row: &SystemMatrixRow, image: &Image) -> Lengthf32 {
    let mut projection = 0.0;
    for (j, w) in system_matrix_row {
        projection += w * image[j]
    }
    projection
}

#[inline]
pub fn back_project(backprojection: &mut [Lengthf32], system_matrix_row: &SystemMatrixRow, projection: Lengthf32) {
    let projection_reciprocal = 1.0 / projection;
    for (j, w) in system_matrix_row {
        backprojection[j] += w * projection_reciprocal;
    }
}
