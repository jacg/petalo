#[derive(Debug, Clone, Copy)]
pub struct Siddon {
    tof: Option<Gaussian>,
}

type Fs<'i> = FoldState<'i, Siddon>;

impl Projector for Siddon {

    fn    project_one_lor<'i>(state: Fs<'i>, lor: &LOR) -> Fs<'i> {
        project_one_lor(state, lor, |projection, lor| ratio_(projection * lor.additive_correction))
    }

    fn sensitivity_one_lor<'i>(state: Fs<'i>, lor: &LOR) -> Fs<'i> {
        project_one_lor(state, lor, |projection, _lor| (-projection).exp())
    }

    fn buffers(fov: FOV) -> SystemMatrixRow {
        let [nx, ny, nz] = fov.n;
        let max_number_of_coupled_voxels_possible = nx + ny + nz - 2;
        SystemMatrixRow(Vec::with_capacity(max_number_of_coupled_voxels_possible))
    }

}

impl Siddon {
    pub fn new(tof: Option<Tof>) -> Self { Self { tof: make_gauss_option(tof) }}
    pub fn notof() -> Self { Self { tof: None } }

    // TODO Should FOV become construction-time argument?
    pub fn new_system_matrix_row(self, lor: &LOR, fov: &FOV) -> SystemMatrixRow {
        let mut system_matrix_row = Self::buffers(*fov);
        match lor_fov_hit(lor, *fov) {
            None => (),
            Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {
                self.update_smatrix_row(
                    &mut system_matrix_row,
                    next_boundary, voxel_size,
                    index, delta_index, remaining,
                    tof_peak,
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
        &self,
        system_matrix_row: &mut SystemMatrixRow,
        mut next_boundary: Vector,
        voxel_size: Vector,
        mut index: i32,
        delta_index: [i32; 3],
        mut remaining: [i32; 3],
        tof_peak: Length,
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
            if let Some(gauss) = &self.tof {
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

fn project_one_lor<'img>(
    state: Fs<'img>,
    lor: &LOR,
    adapt_forward_projection: impl Fn(f32, &LOR) -> f32,
) -> Fs<'img> {
    let Fs { mut backprojection, mut system_matrix_row, image, projector_data } = state;

    // Need to return the state from various match arms
    macro_rules! return_state { () => (return Fs {backprojection, system_matrix_row, image, projector_data}); }

    // Analyse point where LOR hits FOV
    match lor_fov_hit(lor, image.fov) {

        // LOR missed FOV: nothing to be done
        None => return_state!(),

        // Data needed by `system_matrix_elements`
        Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

            // Find non-zero elements (voxels coupled to this LOR) and their values
            Siddon::update_smatrix_row(
                &projector_data,
                &mut system_matrix_row,
                next_boundary, voxel_size,
                index, delta_index, remaining,
                tof_peak,
            );

            // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
            for (i, _) in &system_matrix_row {
                if i >= backprojection.len() { return_state!(); }
            }

            // Sum product of relevant voxels' weights and activities
            let projection = forward_project(&system_matrix_row, image);

            // ... the sum needs to be adapted for the use specific use case:
            // MLEM or sensitivity image generation, are the only ones so far
            let adapted_projection = adapt_forward_projection(projection, lor);

            // Backprojection of LOR onto image
            back_project(&mut backprojection, &system_matrix_row, adapted_projection);
            return_state!();
        }
    }
}

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

// ----- imports ----------------------------------------------------------------------
use units::{
    C, Length, PerLength, Ratio, Time,
    ratio_, mm, mm_,
    in_base_unit,
    uom::ConstZero,
};

use geometry::Vector;

use crate::{
    BoxDim_u, LOR, Point, Pointf32, RatioPoint, RatioVec,
    config::mlem::Tof,
    fov::FOV,
    gauss::{make_gauss_option, Gaussian},
    index::index3_to_1,
    system_matrix::{SystemMatrixRow, forward_project, back_project},
};

use super::{FoldState, Projector};

// ------------------------------ TESTS ------------------------------
#[cfg(test)]
use float_eq::assert_float_eq;

#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use units::{TWOPI, ratio, mm, mm_, uom::ConstZero, Time, todo::Lengthf32};
    use crate::{index::index1_to_3, Point};

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
        let lor = crate::LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0));
        let command = crate::visualize::vislor_command(&fov, &lor);
        println!("\nTo visualize this case, run:\n{}\n", command);

        // Collect voxels traversed by LOR
        let hits = Siddon::notof().new_system_matrix_row(&LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)), &fov);

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

            let summed: Lengthf32 = Siddon::notof().new_system_matrix_row(
                &LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)), &fov
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

#[cfg(test)]
mod sensitivity_image {
    use super::*;
    use units::{mm, mm_, ns, ratio, ratio_, turn, turn_, Angle, Length, Ratio};
    use rstest::{rstest, fixture};
    use float_eq::assert_float_eq;
    use crate::{mlem::mlem, image::Image};

    /// Representation of a straght line in 2D
    /// ax + by + c = 0
    #[derive(Debug, Copy, Clone, PartialEq)]
    struct Line { a: Ratio, b: Ratio, c: Length }

    impl Line {
        /// Construct a line passing through `(x,y)` at `angle` to the positive x-axis.
        fn from_point_and_angle((x,y): (Length, Length), angle: Angle) -> Line {
            let (b,a) = inverse_atan2(angle);
            let b = -b;
            Self { a, b, c: -(a*x + b*y) }
        }

        /// The points at which this line crosses a circle of radius `r`,
        /// centred on the origin.
        fn circle_intersection(self, r: Length) -> Points {
            let eps = mm(000.1) * mm(000.1);
            let Self{ a, b, c } = self;
            let a2_b2 = a*a + b*b;
            let x = - (a*c) / a2_b2;
            let y = - (b*c) / a2_b2;
            if  c*c > r*r*a2_b2        + eps { return Points::Zero }
            if (c*c - r*r*a2_b2).abs() < eps { return Points::One { x,  y } }
            let d = r*r - c*c / a2_b2;
            let m = (d / a2_b2).sqrt();
            let (x1, y1) = (x + m*b, y - m*a);
            let (x2, y2) = (x - m*b, y + m*a);
            Points::Two { x1, y1, x2, y2 }
        }
    }

    impl std::fmt::Display for Line {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let Line { a, b, c } = *self;
            let (a, b, c) = (ratio_(a), ratio_(b), mm_(c));
            write!(f, "{a}x {b:+}y {c:+} = 0")
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq)]
    enum Points { Zero, One{ x: Length, y: Length }, Two{ x1: Length, y1: Length, x2: Length, y2: Length } }

    #[rstest(/**/   x ,   y ,  turns ,      a ,   b ,   c ,
             case( 0.0,  0.0,   0.0  ,     0.0, -1.0,  0.0), // y =  0
             case( 9.0,  0.0,   0.0  ,     0.0, -1.0,  0.0), // y =  0
             case( 0.0,  0.0,   0.25 ,     1.0,  0.0,  0.0), // x =  0
             case( 0.0,  3.0,   0.25 ,     1.0,  0.0,  0.0), // x =  0
             case( 0.0,  2.0,   0.0  ,     0.0, -1.0,  2.0), // y =  2
             case( 5.0,  2.0,   0.0  ,     0.0, -1.0,  2.0), // y =  2
             case( 0.0, -2.0,   0.0  ,     0.0, -1.0, -2.0), // y = -2
             case( 2.0,  0.0,   0.25 ,     1.0,  0.0, -2.0), // x =  2
             case( 2.0,  9.0,   0.25 ,     1.0,  0.0, -2.0), // x =  2
             case(-2.0,  0.0,   0.25 ,     1.0,  0.0,  2.0), // x = -2
             case( 0.0,  0.0,   0.125,     1.0, -1.0,  0.0), // y =  x
             case( 2.0,  0.0,   0.125,     1.0, -1.0, -2.0), // y =  x - 2
             case( 0.0, -2.0,   0.125,     1.0, -1.0, -2.0), // y =  x - 2
             case(-2.0,  0.0,   0.125,     1.0, -1.0,  2.0), // y =  x + 2
             case( 0.0,  2.0,   0.125,     1.0, -1.0,  2.0), // y =  x + 2
             //case( 0.0,  0.0,   0.17618,   1.0, -0.5,  0.0), // y = 2x
    )]
    fn construct_line(x: f32, y: f32, turns: f32, a: f32, b: f32, c: f32) {
        let (x, y, turns) = (mm(x), mm(y), turn(turns));
        let line = Line::from_point_and_angle((x,y), turns);
        let expected = Line { a: ratio(a), b: ratio(b), c: mm(c) };
        println!("constructed: {line}\nexpected   : {expected}");
        assert_float_eq!(ratio_(line.a), a, ulps <= 1);
        assert_float_eq!(ratio_(line.b), b, ulps <= 1);
        assert_float_eq!(   mm_(line.c), c, ulps <= 1);
    }

    #[rstest(/**/   x ,  y , turns,     r   , expected,
             // Vertical lines, close to y = 2
             case( 6.6, 2.1, 0.0  , mm(2.0) , Points::Zero)                                                       ,
             case( 6.6, 2.0, 0.0  , mm(2.0) , Points::One { x : mm( 0.0), y : mm( 2.0)                           }),
             //case( 6.6, 1.9, 0.0  , mm(2.0) , Points::Two { x1: mm(-0.6), y1: mm( 1.9), x2: mm(0.6), y2: mm(1.9) }),
             case( 6.6, 0.0, 0.0  , mm(2.0) , Points::Two { x1: mm(-2.0), y1: mm( 0.0), x2: mm(2.0), y2: mm(0.0) }),
             // Horizontal lines, close to x = 2
             case( 2.1, 9.9, 0.25 , mm(2.0) , Points::Zero)                                                       ,
             case( 2.0, 9.9, 0.25 , mm(2.0) , Points::One { x : mm( 2.0), y : mm( 0.0)                           }),
             //case( 1.9, 9.9, 0.25 , mm(2.0) , Points::Two { x1: mm( 1.9), y1: mm(-0.6), x2: mm(1.9), y2: mm(0.6) }),
             case( 0.0, 9.9, 0.25 , mm(2.0) , Points::Two { x1: mm( 0.0), y1: mm(-2.0), x2: mm(0.0), y2: mm(2.0) }),
    )]
    fn intersect_line_circle(x: f32, y: f32, turns: f32, r: Length, expected: Points) {
        let (x, y, turns) = (mm(x), mm(y), turn(turns));
        let line = Line::from_point_and_angle((x,y), turns);
        let points = line.circle_intersection(r);
        println!("---------> {line} <------------");
        assert_eq!(points, expected);
    }

    /// A version of tan that avoids problems near the poles of tan
    /// angle -> a,b ∈ \[-1, 1\] such that b/a = tan(angle)
    fn inverse_atan2(mut angle: Angle) -> (Ratio, Ratio) {
        let full   : Angle = turn(1.0);
        let half   : Angle = full    / 2.0;
        let quarter: Angle = half    / 2.0;
        let eighth : Angle = quarter / 2.0;
        // Move to preferred range [-1/8, +3/8) turns
        while angle >= 3.0 * eighth {
            angle -= half;
        }
        if angle < eighth { (ratio(1.0)             , angle.tan()) }
        else              { ((quarter - angle).tan(), ratio(1.0) ) }

    }

    /// Tangent of n / 32 full turns
    fn t32(n: i32) -> f32 { ratio_(turn(n as f32 / 32.0).tan()) }

    #[rstest(/**/    turns   ,       x ,       y  ,
             case( 0.0 / 32.0,      1.0,      0.0),
             case( 1.0 / 32.0,      1.0,  t32(1) ),
             case( 2.0 / 32.0,      1.0,  t32(2) ),
             case( 3.0 / 32.0,      1.0,  t32(3) ),
             case( 4.0 / 32.0,      1.0,      1.0),
             case( 5.0 / 32.0,  t32(3) ,      1.0),
             case( 6.0 / 32.0,  t32(2) ,      1.0),
             case( 7.0 / 32.0,  t32(1) ,      1.0),
             case( 8.0 / 32.0,      0.0,      1.0),
             case( 9.0 / 32.0, -t32(1) ,      1.0),
             case(10.0 / 32.0, -t32(2) ,      1.0),
             case(11.0 / 32.0, -t32(3) ,      1.0),
             case(12.0 / 32.0,      1.0,     -1.0),
             case(13.0 / 32.0,      1.0, -t32(3) ),
             case(14.0 / 32.0,      1.0, -t32(2) ),
             case(15.0 / 32.0,      1.0, -t32(1) ),
             case(16.0 / 32.0,      1.0,      0.0),
             case(17.0 / 32.0,      1.0,  t32(1) ),
             case(18.0 / 32.0,      1.0,  t32(2) ),
             case(19.0 / 32.0,      1.0,  t32(3) ),
             case(20.0 / 32.0,      1.0,      1.0),
             case(21.0 / 32.0,  t32(3) ,      1.0),
             case(22.0 / 32.0,  t32(2) ,      1.0),
             case(23.0 / 32.0,  t32(1) ,      1.0),
             case(24.0 / 32.0,      0.0,      1.0),
             case(25.0 / 32.0, -t32(1) ,      1.0),
             case(26.0 / 32.0, -t32(2) ,      1.0),
             case(27.0 / 32.0, -t32(3) ,      1.0),
             case(28.0 / 32.0,      1.0,     -1.0),
             case(29.0 / 32.0,      1.0, -t32(3)),
             case(30.0 / 32.0,      1.0, -t32(2)),
             case(31.0 / 32.0,      1.0, -t32(1)),
    )]
    fn inv_atan2(turns: f32, x: f32, y: f32) {
        let (a,b) = inverse_atan2(turn(turns));
        assert_float_eq!((ratio_(a),ratio_(b)), (x,y), abs <= (3e-7, 4e-7));
    }

    /// `n` uniformly angularly distributed LOR passing through `(x,y)`
    fn n_lors_through(n: usize, (x, y): (Length, Length)) -> Vec<LOR> {
        let mut lors = vec![];
        for angle in n_angles_around_half_circle_starting_from(n, turn(0.01)) {
            let lor = Line::from_point_and_angle((x,y), angle);
            match lor.circle_intersection(DETECTOR_RADIUS) {
                Points::Two { x1, y1, x2, y2 } => {
                    lors.push(LOR::from_components((ns(0.0), ns(0.0)),
                                                   (x1, y1, mm(0.0)),
                                                   (x2, y2, mm(0.0)),
                                                   ratio(1.0)))
                },
                _ => panic!("LOR does not cross detector at two points.")
            }
        }
        lors
    }

    fn n_angles_around_half_circle_starting_from(n: usize, start: Angle) -> impl Iterator<Item = Angle> {
        let mut i = 0;
        let step = turn(0.5) / (n as f32);
        std::iter::from_fn(move || {
            if i < n {
                let angle = start + i as f32 * step;
                i += 1;
                Some(angle)
            } else {
                None
            }
        })
    }

    #[rstest(/**/ n, start    , expected,
             case(0, turn(0.4), vec![]),
             case(1, turn(0.3), vec![turn(0.3)]),
             case(2, turn(0.2), vec![turn(0.2), turn(0.45)]),
             case(2, turn(0.1), vec![turn(0.1), turn(0.35)]),
             case(5, turn(0.2), vec![turn(0.2), turn(0.3), turn(0.4), turn(0.5), turn(0.6)]),
    )]
    fn distributed_angles(n: usize, start: Angle, expected: Vec<Angle>) {
        let angles: Vec<f32> = n_angles_around_half_circle_starting_from(n, start)
            .map(turn_)
            .collect();
        let expected: Vec<f32> = expected
            .into_iter()
            .map(turn_)
            .collect();
        assert_float_eq!(angles, expected, ulps_all <= 1);
    }

    /// Generate points on a square lattice
    fn grid(xs: (i32, i32), ys: (i32, i32)) -> impl Iterator<Item = (i32, i32)> {
        let mut x = xs.0;
        let mut y = ys.0;
        std::iter::from_fn(move || {
            if y > ys.1 { x += 1; y = ys.0 }
            if x > xs.1 { return None }
            y += 1;
            Some((x, (y - 1)))
        })
    }

    // TODO: this should be reused in bin/mlem:main (report_time complicates it)
    /// Return a function which saves images in the given directory
    fn save_each_image_in(directory: String) -> impl FnMut((Image, Osem)) {
        use std::path::PathBuf;
        std::fs::create_dir_all(PathBuf::from(&directory)).unwrap();
        move |(image, Osem{iteration, subset, ..})| {
            let image_path = PathBuf::from(format!("{directory}/{iteration:02}-{subset:02}.raw"));
            crate::io::raw::Image3D::from(&image).write_to_file(&image_path).unwrap();
        }
    }

    /// The position of a rectangular region of interest in the x-y plane, with
    /// an associated count of decays associated with the region
    struct Roi {
        x: (i32, i32),
        y: (i32, i32),
        activity: usize,
    }

    /// Create a vector of coincidences generated by decays in the foreground
    /// and background ROIs. The background only contributes decays in those
    /// voxels which are not covered by any of the foreground ROIs. To
    /// facilitate generating higher statistics, `scale` is applied to each
    /// activity.
    fn trues_from_rois(foreground_rois: &[&Roi], background_roi: &Roi, scale: usize) -> Vec<LOR> {
        let mut activity = std::collections::HashMap::new();
        // Sum activities in voxels covered by foreground_rois
        for roi in foreground_rois {
            for (x,y) in grid(roi.x, roi.y) {
                *activity.entry((x, y)).or_insert(0) += roi.activity * scale;
            }
        }
        // Set activity in background_roi only in voxels which have not been set by the foregrounds
        for (x,y) in grid(background_roi.x, background_roi.y) {
            activity.entry((x, y)).or_insert(background_roi.activity * scale);
        }
        let mut lors = vec![];
        for ((x,y), a) in activity {
            lors.extend(n_lors_through(a, (mm(x as f32), mm(y as f32))))
        }
        lors
    }

    /// Use a rough model to apply a scattering perturbation to `lors`
    ///
    /// The given LORs are adjusted by randomly picking one end and shifting it
    /// around the detector ring by a random angle within `max_deviation`.
    ///
    /// The number of scatters generated per input LOR is controlled by `scale`.
    fn scatter_lors(lors: impl Iterator<Item = LOR>, max_deviation: Angle, scale: i16) -> Vec<LOR> {
        let max_deviation: f32 = turn_(max_deviation);
        use geometry::Point;
        use rand::prelude::*;

        fn rotate_by(delta: Angle, Point { x, y, z }: Point) -> Point {
            let new_angle = y.atan2(x) + delta;
            let x = DETECTOR_RADIUS * new_angle.cos();
            let y = DETECTOR_RADIUS * new_angle.sin();
            Point { x, y, z }
        }

        // How to scatter one LOR
        let scatter = |mut lor: LOR| {
            let mut rng = thread_rng();
            let delta = turn(rng.gen_range(-max_deviation..=max_deviation));
            if rng.gen() { lor.p1 = rotate_by(delta, lor.p1); }
            else         { lor.p2 = rotate_by(delta, lor.p2); }
            lor
        };

        // How to produce multiple scatters derived single LOR
        let scaled_scatter = |lor: LOR| {
            let mut i = 0;
            std::iter::from_fn(move || {
                if i < scale {
                    i += 1;
                    Some(scatter(lor))
                } else {
                    None
                }
            })
        };

        lors.into_iter()
            .flat_map(scaled_scatter)
            .collect()
    }

    use units::in_base_unit;
    pub const DETECTOR_RADIUS: Length = in_base_unit!(50.0);

    // Regions of interest for re-use in tests
    const ACT_1: usize =   0;
    const ACT_2: usize = 100;
    const ACT_3: usize =  80;
    const BG   : usize =  20;
    #[fixture] fn roi_1() -> Roi { Roi { x: (-15,-10), y: (  9,14), activity: ACT_1 } }
    #[fixture] fn roi_2() -> Roi { Roi { x: ( 10, 15), y: (  6,11), activity: ACT_2 } }
    #[fixture] fn roi_3() -> Roi { Roi { x: (- 7,  7), y: ( -8,-4), activity: ACT_3 } }
    #[fixture] fn roi_b() -> Roi { Roi { x: (-20, 20), y: (-20,20), activity: BG    } }

    #[fixture]
    fn fov() -> FOV {
        let n = 51;
        let l = mm(n as f32);
        FOV::new((l, l, mm(1.0)),
                 (n, n,    1   ))
    }

    #[fixture]
    #[once]
    fn trues_and_scatters(roi_1: Roi, roi_2: Roi, roi_3: Roi, roi_b: Roi) -> (Vec<LOR>, Vec<LOR>) {
        // Generate scatters and trues
        let trues = trues_from_rois(&[&roi_1, &roi_2, &roi_3], &roi_b, 1);
        let noise = scatter_lors(trues.iter().cloned(), turn(0.1), 1);

        // // Write LORs to file, for use in testing of CLI scattergram specification
        // let mut hdf5lors: Vec<Hdf5Lor> = vec![];
        // {
        //     let mklor = | &LOR { p1, p2, .. }, prompt | {
        //         let (E1, E2) = match prompt {
        //             Prompt::True    => (511.0, 511.0),
        //             Prompt::Scatter => (450.0, 450.0),
        //             _ => panic!("Not expecting randoms"),
        //         };
        //         Hdf5Lor {
        //             dt: 0.0,
        //             x1: mm_(p1.x), y1: mm_(p1.y), z1: mm_(p1.z),
        //             x2: mm_(p2.x), y2: mm_(p2.y), z2: mm_(p2.z),
        //             q1: f32::NAN, q2: f32::NAN,
        //             E1, E2
        //         }
        //     };
        //     for lor in &trues { hdf5lors.push(mklor(lor, Prompt::True   )); }
        //     for lor in &noise { hdf5lors.push(mklor(lor, Prompt::Scatter)); }
        // }
        // let filename = std::path::PathBuf::from("test-mlem-images/lors.h5");
        // std::fs::create_dir_all(filename.parent().unwrap()).unwrap();
        // hdf5::File::create(filename).unwrap()
        //     .create_group("reco_info").unwrap()
        //     .new_dataset_builder()
        //     .with_data(&hdf5lors)
        //     .create("lors").unwrap();

        (trues, noise)
    }

    use crate::{lorogram::{BuildScattergram as Sc, Prompt}, projector::Siddon, mlem::Osem};

    #[rstest(/**/ name        , correction,
             case("corr-none" , Sc::new()                                        ),
             case("corr-r"    , Sc::new()             .r_bins(20).r_max(mm(30.0))),
             case("corr-phi"  , Sc::new().phi_bins(20)                           ),
             case("corr-r-phi", Sc::new().phi_bins(20).r_bins(20).r_max(mm(30.0))),
    )]
    fn mlem_reco(correction: Sc, name: &str, fov: FOV,
                 trues_and_scatters: &(Vec<LOR>, Vec<LOR>))
    {
        let (trues, noise) = trues_and_scatters;

        let mut sgram = correction.build();

        // Fill scattergam, if required
        if let Some(sgram) = &mut sgram {
            for lor in trues { sgram.fill(Prompt::True   , lor); }
            for lor in noise { sgram.fill(Prompt::Scatter, lor); }
        }

        // Combines trues and scatters into single collection of prompts
        let mut lors = trues.clone();
        lors.extend(noise);

        // Annotate each LOR with additive correction taken from scattergam
        if let Some(sgram) = sgram {
            for mut lor in &mut lors {
                lor.additive_correction = sgram.value(lor);
            }
        }

        // Perform MLEM reconstruction, saving images to disk
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        pool.install(|| {
            mlem(Siddon::notof(), fov, &lors, None, 1)
                .take(10)
                .for_each(save_each_image_in(format!("test-mlem-images/{name}/")));
        });
        //assert!(false);
    }
}
