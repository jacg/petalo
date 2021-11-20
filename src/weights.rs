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

use ncollide3d::query::RayCast;
use ncollide3d::shape::Cuboid;

type Ray      = ncollide3d::query::Ray    <Length>;
type Isometry = ncollide3d::math::Isometry<Length>;
type VecOf<T> = ncollide3d::math::Vector<T>;


// TODO: have another go at getting nalgebra to work with uom.

use crate::types::{BoxDim, Index1, Index3, Index3Weight, Length, Point, Time, Vector, ns_to_mm};
use crate::gauss::make_gauss_option;
use crate::mlem::{index3_to_1, index1_to_3};

const EPS: Length = 1e-5;

// ------------------------------ TESTS ------------------------------
#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use assert_approx_eq::assert_approx_eq;
    use crate::types::TWOPI;

    // --------------------------------------------------------------------------------
    // This set of hand-picked values should be easy to verify by humans. The
    // test performs two checks:
    //
    // 1. The sum of the LOR-lengths within individual voxels equals the
    //    expected total length of LOR in the whole voxel box.
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
    fn hand_picked(p1:   (Length, Length),
                   p2:   (Length, Length),
                   size: (Length, Length),
                   n: (usize, usize),
                   length: Length,
                   expected_voxels: Vec<(usize, usize)>) {

        let p1 = Point::new(p1.0, p1.1, 0.0);
        let p2 = Point::new(p2.0, p2.1, 0.0);
        let vbox = VoxelBox::new((size.0, size.1, 1.0), (n.0, n.1, 1));

        // Values to plug in to visualizer:
        let lor = LOR::new(0.0, 0.0, p1, p2);
        let command = crate::visualize::vislor_command(&vbox, &lor);
        println!("\nTo visualize this case, run:\n{}\n", command);

        // Collect hits
        let hits: Vec<Index3Weight> = LOR::new(0.0, 0.0, p1, p2).active_voxels(&vbox, None, None);

        // Diagnostic output
        for (is, l) in &hits { println!("  ({} {})   {}", is[0], is[1], l) }

        // Check total length through voxel box
        let total_length: Length = hits.iter()
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
    // the total length of the LOR in the voxel box equals the sum of its
    // lengths in the individual voxels.
    proptest! {
        #[test]
        fn sum_of_weights_equals_length_through_box(
            // Activated sensor positions
            r        in  200.0..(300.0 as Length),
            p1_angle in 0.0..(1.0 as Length), // around the circle
            p2_delta in 0.1..(0.9 as Length), // relative to p1_angle
            p1_z     in -200.0..(200.0 as Length),
            p2_z     in -200.0..(200.0 as Length),
            // Voxel box
            dx in  100.0..(150.0 as Length),
            dy in  100.0..(150.0 as Length),
            dz in  100.0..(190.0 as Length),
            nx in  5..50_usize,
            ny in  5..50_usize,
            nz in  5..90_usize,
        ) {
            let p1_theta: Length = p1_angle * TWOPI;
            let p2_theta: Length = p1_theta + (p2_delta * TWOPI);
            let p1 = Point::new(r * p1_theta.cos(), r * p1_theta.sin(), p1_z);
            let p2 = Point::new(r * p2_theta.cos(), r * p2_theta.sin(), p2_z);
            let vbox = VoxelBox::new((dx, dy, dz), (nx, ny, nz));

            // Values to plug in to visualizer:
            let lor = LOR::new(0.0, 0.0, p1, p2);
            let command = crate::visualize::vislor_command(&vbox, &lor);
            println!("\nTo visualize this case, run:\n{}\n", command);

            let summed: Length = LOR::new(0.0, 0.0, p1, p2)
                .active_voxels(&vbox, None, None)
                .iter()
                .inspect(|(i, l)| println!("  ({} {} {}) {}", i[0], i[1], i[2], l))
                .map(|(_index, weight)| weight)
                .sum();

            let a = vbox.entry(&p1, &p2);
            let b = vbox.entry(&p2, &p1);

            let in_one_go = match (a,b) {
                (Some(a), Some(b)) => (a - b).magnitude(),
                _ => 0.0
            };

            assert_approx_eq!(summed, in_one_go, 1e-3);

        }
    }
}

// ---------------------- Implementation -----------------------------------------

// Returned by lor_vbox_hit
//        next_boundary voxel_size index d_index   remaining  tof_peak
type VboxHit = (Vector, Vector,    i32,  [i32; 3], [i32; 3],  Length);

/// Figure out if the LOR hits the voxel box at all. If it does, calculate
/// values needed by find_active_voxels.
#[inline]
pub fn lor_vbox_hit(lor: &LOR, vbox: VoxelBox) -> Option<VboxHit> {

    // Simplify expression of the algorithm by flipping axes so that the
    // direction from p1 to p2 is non-negative along all axes. Remember
    // which directions have been flipped, to recover correct voxel indices.
    let (p1, p2, flipped) = flip_axes(lor.p1, lor.p2);

    // If and where LOR enters voxel box.
    let entry_point: Point = match vbox.entry(&p1, &p2) {
        // If LOR misses the box, immediately return
        None => return None,
        // Otherwise, unwrap the point and continue
        Some(point) => point,
    };

    // How far the entry point is from the TOF peak
    let tof_peak = find_tof_peak(entry_point, p1, p2, lor.dt);

    // Express entry point in voxel coordinates: floor(position) = index of voxel.
    let entry_point = find_entry_point(entry_point, vbox);

    // Bookkeeping information needed during traversal of voxel box
    let (
        index,       // current 1d index into 3d array of voxels
        delta_index, // how the index changes along each dimension
        remaining,   // voxels until edge of vbox in each dimension
    ) = index_trackers(entry_point, flipped, vbox.n);

    // Voxel size in LOR length units: how far must we move along LOR to
    // traverse one voxel, in any dimension.
    let voxel_size = voxel_size(vbox, p1, p2);

    // At what position along LOR is the next voxel boundary, in any dimension.
    let next_boundary = first_boundaries(entry_point, voxel_size);

    // Return the values needed by `find_active_voxels`
    Some((next_boundary, voxel_size, index, delta_index, remaining, tof_peak))
}

/// For a single LOR, place the weights and indices of the active voxels in the
/// `weights` and `indices` parameters. Using output parameters rather than
/// return values, because this function is called in the inner loop, and
/// allocating the vectors of results repeatedly, had a noticeable impact on
/// performance.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn find_active_voxels(
    indices: &mut Vec<usize>,
    weights: &mut Vec<Length>,
    mut next_boundary: Vector,
    voxel_size: Vector,
    mut index: i32,
    delta_index: [i32; 3],
    mut remaining: [i32; 3],
    tof_peak: Length,
    tof: &Option<impl Fn(Length) -> Length>) {

    // How far we have moved since entering the voxel box
    let mut here = 0.0;

    loop {
        // Which voxel boundary will be hit next, and its position
        let (dimension, boundary_position) = next_boundary.argmin();

        // The weight is the length of LOR in this voxel
        let mut weight = boundary_position - here;

        // If TOF enabled, adjust weight
        if let Some(gauss) = &tof {
            weight *= gauss(here - tof_peak);
        }

        // Store the index and weight of the voxel we have just crossed
        if weight > 0.0 {
            indices.push(index as usize);
            weights.push(weight);
        }

        // Move along LOR until it leaves this voxel
        here = boundary_position;

        // Find the next boundary in this dimension
        next_boundary[dimension] += voxel_size[dimension];

        // Move index across the boundary we are crossing
        index += delta_index[dimension];
        remaining[dimension] -= 1;

        // If we have traversed the whole voxel box, we're finished
        if remaining[dimension] == 0 { break; }
    }
}

/// The point at which the LOR enters the voxel box, expressed in a coordinate
/// system with one corner of the box at the origin
#[inline]
fn find_entry_point(mut entry_point: Point, vbox: VoxelBox) -> Vector {
    // Transform coordinates to align box with axes: making the lower boundaries
    // of the box lie on the zero-planes.
    entry_point += vbox.half_width;

    // Express entry point in voxel coordinates: floor(position) = index of voxel.
    let mut entry_point: Vector = entry_point.coords.component_div(&vbox.voxel_size);

    // Floating-point subtractions which should give zero, usually miss very
    // slightly: if this error is negative, the next step (which uses floor)
    // will pick the wrong voxel. Work around this problem by assuming that
    // anything very close to zero is exactly zero.
    entry_point.iter_mut().for_each(|x| if x.abs() < EPS { *x = 0.0 });
    entry_point
}

/// Distance from entry point to the LOR's TOF peak
#[inline]
fn find_tof_peak(entry_point: Point, p1: Point, p2: Point, dt: Time) -> Length {
    let half_lor_length = (p1 - p2).norm() / 2.0;
    let tof_shift = ns_to_mm(dt) / 2.0; // NOTE ignoring refractive index
    let p1_to_peak = half_lor_length - tof_shift;
    let p1_to_entry = (entry_point - p1).norm();
    p1_to_peak - p1_to_entry
}

/// Distances from entry point to the next voxel boundaries, in each dimension
#[inline]
fn first_boundaries(entry_point: Vector, voxel_size: Vector) -> Vector {
    // What fraction of the voxel has already been traversed at the entry
    // point, along any axis.
    let vox_done_fraction: Vector = entry_point - entry_point.map(|x| x.floor());
    // Distances remaining to the nearest boundaries
    (Vector::repeat(1.0) - vox_done_fraction).component_mul(&voxel_size)
}

/// Voxel size expressed in LOR distance units: how far we must move along LOR
/// to cross one voxel in any given dimension. Will be infinite for any axis
/// which is parallel to the LOR.
#[inline]
fn voxel_size(vbox: VoxelBox, p1: Point, p2: Point) -> Vector {
    let lor_direction = (p2-p1).normalize();
    vbox.voxel_size.component_div(&lor_direction)
}

/// Calculate information needed to keep track of progress across voxel box:
/// voxel index and distance remaining until leaving the box
#[inline]
#[allow(clippy::identity_op)]
fn index_trackers(entry_point: Vector, flipped: VecOf<bool>, [nx, ny, nz]: BoxDim) -> (i32, [i32; 3], [i32; 3]) {

    // Find N-dimensional index of voxel at entry point.
    let [ix, iy, iz] = [entry_point.x.floor() as usize,
                        entry_point.y.floor() as usize,
                        entry_point.z.floor() as usize];

    // index is unsigned, but need signed values for delta_index
    let [ix, iy, iz] = [ix as i32, iy as i32, iz as i32];
    let [nx, ny, nz] = [nx as i32, ny as i32, nz as i32];

    // How much the 1d index changes along each dimension
    let delta_index = [
        1       * if flipped[0] { -1 } else { 1 },
        nx      * if flipped[1] { -1 } else { 1 },
        nx * ny * if flipped[2] { -1 } else { 1 },
    ];

    // How many voxels remain before leaving vbox in each dimension
    let remaining = [
        nx - ix ,
        ny - iy ,
        nz - iz ,
    ];

    // 1d index into the 3d arrangement of voxels
    let [ix, iy, iz] = [
        if flipped[0] { nx - 1 - ix } else { ix },
        if flipped[1] { ny - 1 - iy } else { iy },
        if flipped[2] { nz - 1 - iz } else { iz },
    ];
    let index = index3_to_1([ix, iy, iz], [nx, ny, nz]);

    (index, delta_index, remaining)
}

/// Flip axes to ensure that direction from p1 to p2 is non-negative in all
/// dimensions. Return p1 & p2 in flipped coordinate system, along with
/// knowledge of which axes were flipped, so that the indices of subsequently
/// found active voxels can be flipped back into the original coordinate system.
#[inline]
fn flip_axes(mut p1: Point, mut p2: Point) -> (Point, Point, VecOf<bool>) {
    let dimensions = 3;
    let original_lor_direction: Vector = p2 - p1;
    let mut flipped = VecOf::<bool>::repeat(false);
    let mut flip_if_necessary = |n| {
        if original_lor_direction[n] < 0.0 {
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
pub struct VoxelBox {
    pub half_width: Vector,
    pub n: BoxDim,
    pub voxel_size: Vector,
}

impl VoxelBox {

    pub fn new(
        full_size: (Length, Length, Length),
        (nx, ny, nz): (usize, usize, usize)
    ) -> Self {
        let (dx, dy, dz) = full_size;
        let half_width = Vector::new(dx/2.0, dy/2.0, dz/2.0);
        let n = [nx, ny, nz];
        let voxel_size = Self::voxel_size(n, half_width);
            Self { half_width, n, voxel_size, }
    }

    fn voxel_size(n: BoxDim, half_width: Vector) -> Vector {
        // TODO: generalize conversion of VecOf<int> -> VecOf<float>
        let nl: Vector = Vector::new(n[0] as Length, n[1] as Length, n[2] as Length);
        (half_width * 2.0).component_div(&nl)
    }

    pub fn voxel_centre(&self, i: Index3) -> Point {
        //i.map(|n| n as f64 + 0.5).component_mul(&self.voxel_size).into()
        let s = self.voxel_size;
        Point::new((i[0] as Length + 0.5) * s.x - self.half_width[0],
                   (i[1] as Length + 0.5) * s.y - self.half_width[1],
                   (i[2] as Length + 0.5) * s.z - self.half_width[2],)
    }

    pub fn voxel_centre1(&self, i: Index1) -> Point {
        self.voxel_centre(index1_to_3(i, self.n))
    }

    pub fn entry(&self, p1: &Point, p2: &Point) -> Option<Point> {
        let lor_direction: Vector = (p2 - p1).normalize();
        let lor_length   : Length = (p2 - p1).norm();
        let lor: Ray = Ray::new(*p1, lor_direction);
        let iso: Isometry = Isometry::identity();
        Cuboid::new(self.half_width)
            .toi_with_ray(&iso, &lor, lor_length, true)
            .map(|toi| lor.origin + lor.dir * toi)
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
    fn test_voxel_centre(index: Index3, expected_position: [Length; 3]) {
        let vbox = VoxelBox::new((4.0, 4.0, 4.0), (2,2,2));
        let c = vbox.voxel_centre(index);
        assert!([c.x, c.y, c.z] == expected_position);
    }
}

//--------------------------------------------------------------------------------
/// Line Of Response: 2 spacetime vectors indicating the positions and times of
/// coincident detector element activations
#[derive(Clone, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct LOR {
    pub p1: Point,
    pub p2: Point,
    pub dt: Time,
}

impl LOR {
    pub fn new(t1: Time, t2: Time, p1: Point, p2: Point) -> Self {
        Self { p1, p2, dt: t2 - t1 }
    }

    pub fn from_components(t1: Time, t2: Time,
                           x1: Length, y1: Length, z1: Length,
                           x2: Length, y2: Length, z2: Length) -> Self
    {
        Self::new(t1, t2, Point::new(x1,y1,z1), Point::new(x2,y2,z2))
    }

    pub fn active_voxels(&self, vbox: &VoxelBox, cutoff: Option<Length>, sigma: Option<Length>) -> Vec<Index3Weight> {

        let tof = make_gauss_option(sigma, cutoff);
        let mut weights = vec![];
        let mut indices = vec![];
        match lor_vbox_hit(self, *vbox) {
            None => (),
            Some((next_boundary, voxel_size, index, delta_index, remaining, tof_peak)) => {
                find_active_voxels(
                    &mut indices, &mut weights,
                    next_boundary, voxel_size,
                    index, delta_index, remaining,
                    tof_peak, &tof
                );

            }
        }
        indices.into_iter()
            .map(|i| index1_to_3(i, vbox.n))
            .zip(weights.into_iter())
            .collect()
    }
}

use core::fmt;
impl fmt::Display for LOR {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (p, q) = (self.p1, self.p2);
        write!(f, "<LOR ({:8.2} {:8.2} {:8.2}) ({:8.2} {:8.2} {:8.2}) {:7.2}ns {:7.2}mm /{:7.2} >",
               p.x, p.y, p.z,
               q.x, q.y, q.z,
               self.dt, ns_to_mm(self.dt) / 2.0,
               (p-q).norm()
        )
    }
}
