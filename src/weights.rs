use ncollide3d as nc;
use nc::query::RayCast;

use ndarray::azip;

// TODO: have another go at getting nalgebra to work with uom.

#[allow(clippy::excessive_precision)] // Precision needed when features=f64
pub const C: Length = 0.299792458; // mm / ps
pub const DISTANCE_AS_TOF_DELTA: Length = 2.0 / C;
pub const TOF_DT_AS_DISTANCE: Ratio = C / 2.0;

// TODO: no thought has been given to what should be public or private.

#[cfg    (feature = "f64") ] pub const PRECISION: u8 = 64;
#[cfg(not(feature = "f64"))] pub const PRECISION: u8 = 32;

#[cfg    (feature = "f64") ] pub type Length = f64;
#[cfg(not(feature = "f64"))] pub type Length = f32;

#[cfg    (feature = "f64") ] const EPS: Length = 1e-14;
#[cfg(not(feature = "f64"))] const EPS: Length = 1e-5;

pub type Time   = Length;
pub type Weight = Length;
pub type Ratio = Length;

type Vector = nc::math ::Vector<Length>;
pub type Point  = nc::math ::Point <Length>;
type Ray    = nc::query::Ray   <Length>;
type Isometry = nc::math::Isometry<Length>;

type VecOf<T> = nc::math::Vector<T>;

pub type Index3 = [usize; 3];
pub type Index1 = usize;
type BoxDim = [usize; 3];

//type Index1Weight = (Index1, Weight);
type Index3Weight = (Index3, Weight);

const TWOPI: Length = std::f64::consts::TAU as Length;

// This algorithm is centred around two key simplifications:
//
// 1. Express the voxel size in terms of the LOR components. This allows trivial
//    calculation of how far we must move along the LOR before reaching a voxel
//    boundary, in any dimension.
//
// 2. Exploit symmetry to simplify dealing with directions: flip axes so that
//    the direction of the LOR has non-negative components. The algorithm can
//    then assume that all progress is in the positive direction. Any voxels
//    indices calculated by the algorithm, must be flipped back to the original
//    coordinate system.

// ------------------------------ TESTS ------------------------------
#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use assert_approx_eq::assert_approx_eq;

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

            #[cfg    (feature = "f64") ] assert_approx_eq!(summed, in_one_go);
            #[cfg(not(feature = "f64"))] assert_approx_eq!(summed, in_one_go, 1e-3);

        }
    }
}
//--------------------------------------------------------------------------------
type Intensity = Length;

type ImageData = Vec<Intensity>;

#[derive(Clone)]
pub struct Image {
    vbox: VoxelBox,
    pub data: ImageData,
    size: usize,
}

impl core::ops::IndexMut<Index1> for Image {
    fn index_mut(&mut self, i: Index1) -> &mut Self::Output { &mut self.data[i] }
}

impl core::ops::Index<Index1> for Image {
    type Output = Intensity;
    fn index(&self, i: Index1) -> &Self::Output { &self.data[i] }
}

impl Image {

    pub fn mlem<'a>(vbox: VoxelBox,
                    measured_lors: &'a [LOR],
                    sigma        :     Option<Time>,
                    cutoff       :     Option<Ratio>,
                    smatrix      : &'a ImageData,
    ) -> impl Iterator<Item = Image> + 'a {

        // Start off with a uniform image
        let mut image = Self::ones(vbox);

        // Return an iterator which generates an infinite sequence of images,
        // each one made by performing one MLEM iteration on the previous one
        std::iter::from_fn(move || {
            image.one_iteration(measured_lors, smatrix, sigma, cutoff);
            Some(image.clone()) // TODO see if we can sensibly avoid cloning
        })
    }

    fn one_iteration(&mut self, measured_lors: &[LOR], smatrix: &ImageData, sigma: Option<Time>, cutoff: Option<Ratio>) {

        // Accumulator for all backprojection contributions in this iteration
        let mut backprojection = self.zeros_buffer();

        // TOF adjustment to apply to the weights
        let tof: Option<_> = make_gauss_option(sigma, cutoff);

        // Storage space for the weights and indices of the active voxels
        // (allocating new result vectors for each LOR had a noticeable runtime cost)
        let (mut weights, mut indices) = {
            let [nx, ny, nz] = self.vbox.n;
            let max_number_of_active_voxels_possible = nx + ny + nz - 2;
            (Vec::with_capacity(max_number_of_active_voxels_possible),
             Vec::with_capacity(max_number_of_active_voxels_possible))
        };

        // For each measured LOR ...
        measured_lors.iter().for_each(|lor| {

            // Analyse point where LOR hits voxel box
            match lor_vbox_hit(lor, self.vbox) {

                // LOR missed voxel box: nothing to be done
                None => return,

                // Data needed by `find_active_voxels`
                Some((next_boundary, voxel_size, index, delta_index, remaining, tof_peak)) => {

                    // Throw away previous search results
                    weights.clear();
                    indices.clear();

                    // Find active voxels and their weights
                    find_active_voxels(
                        &mut indices, &mut weights,
                        next_boundary, voxel_size,
                        index, delta_index, remaining,
                        tof_peak, &tof
                    );

                    // Forward projection of current image into this LOR
                    let projection = forward_project(&weights, &indices, self);

                    // Backprojection of LOR onto image
                    back_project(&mut backprojection, &weights, &indices, projection);
                }
            }
        });

        apply_sensitivity_matrix(&mut self.data, &backprojection, smatrix);

    }

    // A new empty data store with matching size
    fn zeros_buffer(&self) -> ImageData { vec![0.0; self.size] }
    pub fn ones(vbox: VoxelBox) -> Self {
        let [x,y,z] = vbox.n;
        let size = x * y * z;
        Self { data: vec![1.0; size], vbox, size}
    }

}

fn lor_vbox_hit(lor: &LOR, vbox: VoxelBox)
                -> Option<(Vector, Vector, i32, [i32; 3], [i32; 3], Length)>
{

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
    let tof_peak = find_tof_peak(entry_point, p1, p2, lor.t1, lor.t2);

    // Express entry point in voxel coordinates: floor(position) = index of voxel.
    let entry_point = align_entry_point(entry_point, vbox);

    // Bookkeeping information needed during traversal of voxel box
    let (
        index,       // current 1d index into 3d array of voxels
        delta_index, // how the index changes along each dimension
        remaining,   // voxels until edge of vbox in each dimension
    ) = index_trackers(entry_point, flipped, vbox.n);

    // Voxel size in LOR length units: how far must we move along LOR to
    // traverse one voxel, in any dimension.
    let voxel_size = voxel_size(vbox, p1, p2);

    // How far we must travel along LOR before hitting next voxel boundary,
    // in any dimension.
    let next_boundary = first_boundaries(entry_point, voxel_size);

    // Return the values needed by `find_active_voxels`
    Some((next_boundary, voxel_size, index, delta_index, remaining, tof_peak))
}

fn find_active_voxels(
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

        // Move along LOR until it leaves this voxel
        here = boundary_position;

        // Find the next boundary in this dimension
        next_boundary[dimension] += voxel_size[dimension];

        // Move index across the boundary we are crossing
        let previous_index = index;
        index += delta_index[dimension];
        remaining[dimension] -= 1;

        // If we have traversed the whole voxel box
        if remaining[dimension] == 0 { break; }

        // Store the index of the voxel we have just crossed, along with
        // the distance that the LOR covered in that voxel.
        if weight > 0.0 {
            indices.push(previous_index as usize);
            weights.push(weight);
        }
    }
}

#[inline]
fn align_entry_point(mut entry_point: Point, vbox: VoxelBox) -> Vector {
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

#[inline]
fn find_tof_peak(entry_point: Point, p1: Point, p2: Point, t1: Time, t2: Time) -> Length {
    let p1_to_peak = 0.5 * ((p1 - p2).norm() + C * (t1 - t2));
    let p1_to_entry = (entry_point - p1).norm();
    p1_to_peak - p1_to_entry
}

#[inline]
fn first_boundaries(entry_point: Vector, voxel_size: Vector) -> Vector {
    // What fraction of the voxel has already been traversed at the entry
    // point, along any axis.
    let vox_done_fraction: Vector = entry_point - entry_point.map(|x| x.floor());
    // Distances remaining to the nearest boundaries
    (Vector::repeat(1.0) - vox_done_fraction).component_mul(&voxel_size)
}

#[inline]
fn voxel_size(vbox: VoxelBox, p1: Point, p2: Point) -> Vector {
    let lor_direction = (p2-p1).normalize();
    vbox.voxel_size.component_div(&lor_direction)
}

#[inline]
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

#[inline]
fn forward_project(weights: &Vec<Length>, indices: &Vec<usize>, image: &Image) -> Length {
    let mut projection = 0.0;
    for (w, j) in weights.iter().zip(indices.iter()) {
        projection += w * image[*j]
    }
    projection
}

#[inline]
fn back_project(backprojection: &mut Vec<Length>, weights: &Vec<Length>, indices: &Vec<usize>, projection: Length) {
    let projection_reciprocal = 1.0 / projection;
    for (w, j) in weights.iter().zip(indices.iter()) {
        backprojection[*j] += w * projection_reciprocal;
    }
}


fn apply_sensitivity_matrix(image: &mut ImageData, backprojection: &Vec<Length>, smatrix: &ImageData) {
    //  TODO express with Option<matrix> and mul reciprocal
    // Apply Sensitivity matrix
    azip!((voxel in image, &b in backprojection, &s in smatrix) {
        if s > 0.0 { *voxel *= b / s }
        else       { *voxel  = 0.0   }
    })
}

//--------------------------------------------------------------------------------
#[derive(Clone, Copy, Debug)]
pub struct VoxelBox {
    pub half_width: Vector,
    pub n: BoxDim,
    pub voxel_size: Vector,
}

impl VoxelBox {

    pub fn new(
        (dx, dy, dz): (Length, Length, Length),
        (nx, ny, nz): (usize, usize, usize)
    ) -> Self {
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
        Point::new((i[0] as Length + 0.5) * s.x,
                   (i[1] as Length + 0.5) * s.y,
                   (i[2] as Length + 0.5) * s.z,)
    }

    pub fn entry(&self, p1: &Point, p2: &Point) -> Option<Point> {
        let lor_direction: Vector = (p2 - p1).normalize();
        let lor_length   : Length = (p2 - p1).norm();
        let lor: Ray = Ray::new(*p1, lor_direction);
        let iso: Isometry = Isometry::identity();
        nc::shape::Cuboid::new(self.half_width)
            .toi_with_ray(&iso, &lor, lor_length, true)
            .map(|toi| lor.origin + lor.dir * toi)
    }

}

//--------------------------------------------------------------------------------

fn make_gauss(sigma: Length, cutoff: Option<Length>) -> impl Fn(Length) -> Length {
    let root_two_pi = TWOPI.sqrt() as Length;
    let peak_height = 1.0 / (sigma * root_two_pi);
    let cutoff = cutoff.map_or(std::f32::INFINITY as Length, |width| width * sigma);
    move |x| {
        if x.abs() < cutoff {
            let y = x / sigma;
            let z = y * y;
            peak_height * (-0.5 * z).exp()
        } else {
            0.0
        }
    }
}

pub fn make_gauss_option(sigma: Option<Length>, cutoff: Option<Length>) -> Option<impl Fn(Length) -> Length> {
    sigma.map(|sigma| make_gauss(sigma * TOF_DT_AS_DISTANCE, cutoff))
}

//--------------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct LOR {
    pub p1: Point,
    pub p2: Point,
    pub t1: Time,
    pub t2: Time,
}

impl LOR {
    pub fn new(t1: Time, t2: Time, p1: Point, p2: Point) -> Self { Self { t1, t2, p1, p2 } }

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


// --------------------------------------------------------------------------------
//                  Conversion between 1d and 3d indices

use std::ops::{Add, Div, Mul , Rem};

fn index3_to_1<T>([ix, iy, iz]: [T; 3], [nx, ny, _nz]: [T; 3]) -> T
where
    T: Mul<Output = T> + Add<Output = T>
{
    ix + (iy + iz * ny) * nx
}

fn index1_to_3<T>(i: T, [nx, ny, _nz]: [T; 3]) -> [T; 3]
where
    T: Mul<Output = T> +
    Div<Output = T> +
    Rem<Output = T> +
    Copy
{
    let z = i / (nx * ny);
    let r = i % (nx * ny);
    let y = r / nx;
    let x = r % nx;
    [x,y,z]
}


#[cfg(test)]
mod test_index_conversion {
    use super::*;
    use rstest::rstest;

    // -------------------- Some hand-picked examples ------------------------------
    #[rstest(/**/    size   , index3 , index1,
             // 1-d examples
             case([ 1, 1, 1], [0,0,0],   0),
             case([ 9, 1, 1], [3,0,0],   3),
             case([ 1, 8, 1], [0,4,0],   4),
             case([ 1, 1, 7], [0,0,5],   5),
             // Counting in binary: note digit reversal
             case([ 2, 2, 2], [0,0,0],   0),
             case([ 2, 2, 2], [1,0,0],   1),
             case([ 2, 2, 2], [0,1,0],   2),
             case([ 2, 2, 2], [1,1,0],   3),
             case([ 2, 2, 2], [0,0,1],   4),
             case([ 2, 2, 2], [1,0,1],   5),
             case([ 2, 2, 2], [0,1,1],   6),
             case([ 2, 2, 2], [1,1,1],   7),
             // Relation to decimal: note reversal
             case([10,10,10], [1,2,3], 321),
             case([10,10,10], [7,9,6], 697),
    )]
    fn hand_picked(size: Index3, index3: Index3, index1: usize) {
        assert_eq!(index3_to_1(index3, size), index1);
        assert_eq!(index1_to_3(index1, size), index3);
    }

    // -------------------- Exhaustive roundtrip testing ------------------------------
    use proptest::prelude::*;

    // A strategy that picks 3-d index limits, and a 1-d index guaranteed to lie
    // within those bounds.
    fn size_and_in_range_index() -> impl Strategy<Value = (Index3, usize)> {
        [1..200_usize, 1..200_usize, 1..200_usize]
            .prop_flat_map(|i| (Just(i), 1..(i[0] * i[1] * i[2])))
    }

    proptest! {
        #[test]
        fn index_roundtrip((size, index) in size_and_in_range_index()) {
            let there = index1_to_3(index, size);
            let back  = index3_to_1(there, size);
            assert_eq!(back, index)
        }

    }
}
