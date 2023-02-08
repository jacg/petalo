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

use units::{C, Length, Ratio, PerLength, Time, todo::Weightf32, in_base_unit};
use units::{mm, mm_, ns_, ratio_};
use units::todo::Lengthf32;
use crate::Index3_u;
use crate::{Point, Vector, RatioPoint, RatioVec};
use crate::fov::FOV;

use crate::gauss::make_gauss_option;
use crate::index::index1_to_3;
use crate::config::mlem::Tof;

// ------------------------------ TESTS ------------------------------
#[cfg(test)]
use float_eq::assert_float_eq;

#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use units::{TWOPI, ratio};

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
        //let hits: SystemMatrixRow = LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)).active_voxels(&fov, None);
        let hits = SystemMatrixRow::new(&LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)), &fov, None);

        // Diagnostic output
        for (is, l) in &hits.0 { println!("  ({} {})   {}", is[0], is[1], l) }

        // Check total length through FOV
        let total_length: Lengthf32 = hits.0.iter()
            .map(|(_index, weight)| weight)
            .sum();
        assert_float_eq!(total_length, length, ulps <= 1);

        // Check voxels hit
        let voxels: Vec<(usize, usize)> = hits.0.into_iter()
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
            let (r, p1_z, p2_z) = (mm(r), mm(p1_z), mm(p2_z));
            let p1_theta = p1_angle * TWOPI;
            let p2_theta = p1_theta + (p2_delta * TWOPI);
            let p1 = Point::new(r * p1_theta.cos(), r * p1_theta.sin(), p1_z);
            let p2 = Point::new(r * p2_theta.cos(), r * p2_theta.sin(), p2_z);
            let fov = FOV::new((mm(dx), mm(dy), mm(dz)), (nx, ny, nz));

            // Values to plug in to visualizer:
            let lor = LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0));
            let command = crate::visualize::vislor_command(&fov, &lor);
            println!("\nTo visualize this case, run:\n{}\n", command);

            let summed: Lengthf32 = SystemMatrixRow::new(
                &LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)), &fov, None
            ).into_iter()
             .inspect(|(i, l)| println!("  ({} {} {}) {}", i[0], i[1], i[2], l))
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
pub type SystemMatrixElement = (Index3_u, Weightf32);

pub struct SystemMatrixRow(Vec<SystemMatrixElement>);

impl SystemMatrixRow {
    pub fn new(lor: &LOR, fov: &FOV, tof: Option<Tof>) -> Self {
        use crate::fov::{lor_fov_hit, FovHit};
        let tof = make_gauss_option(tof);
        let mut weights = vec![];
        let mut indices = vec![];
        match lor_fov_hit(lor, *fov) {
            None => (),
            Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {
                system_matrix_elements(
                    &mut indices, &mut weights,
                    next_boundary, voxel_size,
                    index, delta_index, remaining,
                    tof_peak, &tof
                );

            }
        }
        SystemMatrixRow(indices.into_iter()
            .map(|i| index1_to_3(i, fov.n))
            .zip(weights.into_iter())
            .collect())
    }
    pub fn iter(&self) -> std::slice::Iter<SystemMatrixElement> { self.0.iter() }
}

impl IntoIterator for SystemMatrixRow {
    type Item = SystemMatrixElement;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
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
    mut next_boundary: Vector,
    voxel_size: Vector,
    mut index: i32,
    delta_index: [i32; 3],
    mut remaining: [i32; 3],
    tof_peak: Length,
    tof: &Option<impl Fn(Length) -> PerLength>) {

    // How far we have moved since entering the FOV
    let mut here = Length::ZERO;

    loop {
        // Which voxel boundary will be hit next, and its position
        let (dimension, boundary_position) = next_boundary.argmin();

        // The weight is the length of LOR in this voxel
        let mut weight = boundary_position - here;

        // If TOF enabled, adjust weight
        if let Some(gauss) = &tof {
            let g: PerLength = gauss(here - tof_peak);
            // TODO Normalization
            let completely_arbitrary_factor = 666.0;
            let g: f32 = ratio_(mm(completely_arbitrary_factor) * g);
            weight *= g;
        }

        // Store the index and weight of the voxel we have just crossed
        if weight > Length::ZERO {
            indices.push(index as usize);
            weights.push(mm_(weight));
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

use units::uom::ConstZero;

const EPS: Ratio = in_base_unit!(1e-5);

/// The point at which the LOR enters the FOV, expressed in a coordinate
/// system with one corner of the FOV at the origin.
#[inline]
pub fn find_entry_point(entry_point: Point, fov: FOV) -> RatioPoint {
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
pub fn find_tof_peak(entry_point: Point, p1: Point, p2: Point, dt: Time) -> Length {
    let half_lor_length = (p1 - p2).norm() / 2.0;
    let tof_shift = C * dt / 2.0; // NOTE ignoring refractive index
    let p1_to_peak = half_lor_length - tof_shift;
    let p1_to_entry = (entry_point - p1).norm();
    p1_to_peak - p1_to_entry
}


/// Distances from entry point to the next voxel boundaries, in each dimension
#[inline]
pub fn first_boundaries(entry_point: RatioPoint, voxel_size: Vector) -> Vector {
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
pub fn voxel_size(fov: FOV, p1: Point, p2: Point) -> Vector {
    // TODO: The units are a bit dodgy here. See the TODOs for
    // Vector::{normalize,component_div}
    let lor_direction = (p2-p1).normalize();
    fov.voxel_size.component_div(lor_direction)
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

}

use core::fmt;
impl fmt::Display for LOR {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (p, q) = (self.p1, self.p2);
        write!(f, "<LOR ({:8.2} {:8.2} {:8.2}) ({:8.2} {:8.2} {:8.2}) {:7.2}ns {:7.2}mm /{:7.2} >",
               mm_(p.x), mm_(p.y), mm_(p.z),
               mm_(q.x), mm_(q.y), mm_(q.z),
               ns_(self.dt), mm_(self.dt * C) / 2.0,
               mm_((p-q).norm())
        )
    }
}
