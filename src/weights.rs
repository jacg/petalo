use ncollide3d as nc;
use nc::query::RayCast;

use ndarray as nda;

// TODO: have another go at getting nalgebra to work with uom.
const c : Length = 3e3; // cm / ns

// TODO: no thought has been given to what should be public or private.


pub type Length = f64;
pub type Time = f64;
pub type Weight = f64;

type Vector = nc::math ::Vector<Length>;
pub type Point  = nc::math ::Point <Length>;
type Ray    = nc::query::Ray   <Length>;
type Isometry = nc::math::Isometry<Length>;

type VecOf<T> = nc::math::Vector<T>;

pub type Index3 = [usize; 3];
pub type Index1 = usize;
type BoxDim = VecOf<usize>;

// TODO: propably want to use Index1 as much as possible
type Index3Weight = (Index3, Weight);
type Index1Weight = (Index1, Weight);

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

/// An iterator which yields the N-dimensional indices of voxels which have been
/// traversed by a LOR, along with the distance the LOR covered inside that
/// voxel.
pub enum WeightsAlongLOR {

    // We are at a point outside of the voxel box: no bookkeeping to be done.
    Outside,

    // We are traversing the voxel box:
    Inside {

        // How many LOR distance units must be travelled before reaching the
        // next voxel boundary in any dimension.
        to_boundary: Vector,

        // Voxel box size, as number of voxels in each dimension
        n_voxels: BoxDim,

        // Dimensions of the voxels expressed in LOR distance units. Used to
        // reset components of `to_boundary` when they reach 0.
        voxel_size: Vector,

        // The flipped index of the voxel we have just entered. Must be flipped
        // back before yielding to client.
        index: Index3,

        // We exploit the symmetries of the system by flipping some axes. Here
        // we record which axes have been flipped and must be adjusted when
        // calculating indices to be yielded.
        flipped: VecOf<bool>,
    }
}

impl WeightsAlongLOR {
    pub fn new(mut p1: Point, mut p2: Point, vbox: &VoxelBox) -> Self {

        // This function works in an arbitrary number of dimensions. In order to
        // iterate over all dimensions, we need to know how many there are.
        let dimensions = p1.len();

        // Simplify expression of the algorithm by flipping axes so that the
        // direction from p1 to p2 is non-negative along all axes. Remember
        // which directions have been flipped, to recover correct voxel indices.
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

        // Find if and where LOR enters voxel box.
        let mut entry_point: Point = match vbox.entry(&p1, &p2) {
            // If LOR misses the box, immediately return an iterator which will
            // generate no hits.
            None => return Self::Outside,
            // Otherwise, unwrap the point and continue calculating a more
            // detailed iterator state
            Some(point) => point,
        };

        // Transform coordinates to align box with axes: making the lower
        // boundaries of the box lie on the zero-planes.
        entry_point += vbox.half_width;

        // Express entry point in voxel coordinates: floor(position) = index of voxel.
        let mut entry_point: Vector = entry_point.coords.component_div(&vbox.voxel_size);

        // Floating-point subtractions which should give zero, usually miss very
        // slightly: if this error is negative, the next step (which uses floor)
        // will pick the wrong voxel. Work around this problem by assuming that
        // anything very close to zero is exactly zero.
        entry_point.iter_mut().for_each(|x| if x.abs() < 1e-7 { *x = 0.0 });

        // Find N-dimensional index of voxel at entry point.
        let index: Index3 = [entry_point.x.floor() as usize,
                             entry_point.y.floor() as usize,
                             entry_point.z.floor() as usize];//entry_point.map(|x| x.floor() as usize);

        // Voxel size in LOR length units: how far must we move along LOR to
        // traverse one voxel, in any dimension.
        let lor_direction = (p2-p1).normalize();
        let voxel_size: Vector = vbox.voxel_size.component_div(&lor_direction);

        // What fraction of the voxel has already been traversed at the entry
        // point, along any axis.
        let vox_done_fraction: Vector = entry_point - entry_point.map(|x| x.floor());

        // How far we must travel along LOR before hitting next voxel boundary,
        // in any dimension.
        let to_boundary: Vector = (Vector::repeat(1.0) - vox_done_fraction)
            .component_mul(&voxel_size);

        // Initial iterator state: the point where LOR enters voxel box (in
        // voxel coordinates), along with bookkeeping information.
        Self::Inside { to_boundary, voxel_size, n_voxels: vbox.n, flipped, index }
    }
}

impl Iterator for WeightsAlongLOR {

    // Generate one item for each voxel crossed by the LOR. Each item contains
    // the N-dimensional index of the voxel, and the length of the LOR within
    // that voxel.
    type Item = Index3Weight;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            // If we are not inside the voxel box, then either the LOR
            // completely missed the box, or we have traversed it and come out
            // of the other side. In either case, there is nothing more to be
            // done: The iteration is finished.
            Self::Outside => None,

            // If we are inside the voxel box, then we are at the boundary of a voxel
            Self::Inside { flipped, index, n_voxels, to_boundary, voxel_size } => {

                // Remember index of the voxel we are about to cross (flipped
                // back from our algorithm's internal coordinate system, to the
                // client's original coordinate system).
                let mut true_index = [0; 3];
                for n in 0..3 {
                    if flipped[n] { true_index[n] = n_voxels[n] - 1 - index[n]; }
                    else          { true_index[n] =                   index[n]; }
                }

                // Which boundary will be hit next, and how soon
                let (dimension, distance) = to_boundary.argmin();

                // Move along LOR until we hit voxel boundary
                *to_boundary -= Vector::repeat(distance);

                // For any dimension in which we have reached a voxel boundary
                for dimension in 0..index.len() {
                    if to_boundary[dimension] <= 0.0 {
                        // Reset distance to next boundary
                        to_boundary[dimension] = voxel_size[dimension];
                        // Move index into next voxel
                        index[dimension] += 1;
                    }
                }

                // If we have traversed the whole voxel box
                if index[dimension] >= n_voxels[dimension] {
                    // no more steps need to be taken after this one.
                    *self = Self::Outside
                }

                // Yield the N-dimensional index of the voxel we have just
                // crossed (expressed in the client's coordinate system), along
                // with the distance that the LOR covered in that voxel.
                Some((true_index, distance))
            }
        }
    }
    // TODO: iterator hints
}

// ------------------------------ TESTS ------------------------------
#[cfg(test)]
mod test {
    use super::*;
    #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};
    use rstest::rstest;
    use assert_approx_eq::assert_approx_eq;

    const TWOPI: Length = std::f64::consts::TAU as Length;

    // --------------------------------------------------------------------------------
    // This set of hand-picked values should be easy to verify by humans. The
    // test performs two checks:
    //
    // 1. The sum of the LOR-lengths within individual voxels equals the
    //    expected total length of LOR in the whole voxel box.
    //
    // 2. The indices of the voxels traversed by the LOR are as expected.
    #[rstest(/**/         p1        ,    p2       ,    size   ,   n  ,  length  , expected_voxels,
             // symmetric 3x3, diagonal LOR under all four axis flip combinations
             case((-30.0, -30.0), ( 30.0, 30.0), (5.0, 5.0), (3,3), 14.142135, vec![(0,0), (1,1), (2,2)]),
             case(( 30.0, -30.0), (-30.0, 30.0), (5.0, 5.0), (3,3), 14.142135, vec![(2,0), (1,1), (0,2)]),
             case((-30.0,  30.0), ( 30.0,-30.0), (5.0, 5.0), (3,3), 14.142135, vec![(0,2), (1,1), (2,0)]),
             case(( 30.0,  30.0), (-30.0,-30.0), (5.0, 5.0), (3,3), 14.142135, vec![(2,2), (1,1), (0,0)]),
             // like case 1, but with asymmetric voxels
             case((-30.0, -30.0), ( 30.0, 30.0), (5.0, 5.0), (3,2), 14.142135, vec![(0,0), (1,0), (1,1), (2,1)]),
             case((-30.0, -30.0), ( 30.0, 30.0), (5.0, 5.0), (2,3), 14.142135, vec![(0,0), (0,1), (1,1), (1,2)]),
             // vertical / horizontal off-centre LOR
             case((  5.4, -20.0), (  5.4, 10.0), (5.5, 4.5), (9,4),  9.0     , vec![(8,0), (8,1), (8,2), (8,3)]),
             case((-15.0,  -4.0), ( 15.0, -4.0), (4.0, 5.0), (4,3),  8.0     , vec![(0,0), (1,0), (2,0), (3,0)]),
    )]
    fn hand_picked(p1:   (Length, Length),
                   p2:   (Length, Length),
                   size: (Length, Length),
                   n: (usize, usize),
                   length: Length,
                   expected_voxels: Vec<(usize, usize)>) {

        let p1 = Point::new(p1.0, p1.1, 0.0);
        let p2 = Point::new(p2.0, p2.1, 0.0);
        let vbox = VoxelBox::new((size.0, size.1, 0.0), (n.0, n.1, 1));

        //crate::visualize::lor_weights(p1, p2, vbox.clone());

        // Collect hits
        let hits: Vec<Index3Weight> =
            WeightsAlongLOR::new(p1, p2, &vbox)
            .inspect(|(is, l)| println!("  ({} {})   {}", is[0], is[1], l))
            .collect();

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

            // // Values to plug in to visualizer:
            // println!("let p1 = Point::new({}, {}, {});", p1.x, p1.y, p1.z);
            // println!("let p2 = Point::new({}, {}, {});", p2.x, p2.y, p2.z);
            // println!("let vbox = VoxelBox::new(({}, {}, {}), ({}, {}, {}));",
            //          vbox.aabb.half_extents.x, vbox.aabb.half_extents.y, vbox.aabb.half_extents.z,
            //          vbox.n.x, vbox.n.y, vbox.n.z);

            let summed: Length = WeightsAlongLOR::new(p1, p2, &vbox)
                .inspect(|(i, l)| println!("  ({} {} {}) {}", i[0], i[1], i[2], l))
                .map(|(_index, weight)| weight)
                .sum();

            let a = vbox.entry(&p1, &p2);
            let b = vbox.entry(&p2, &p1);

            let in_one_go = match (a,b) {
                (Some(a), Some(b)) => (a - b).magnitude(),
                _ => 0.0
            };

            assert_approx_eq!(summed, in_one_go);
        }
    }
}
//--------------------------------------------------------------------------------
type Intensity = f64;

//type ImageData = nda::Array3<Intensity>;
type ImageData = Vec<Intensity>;

pub struct Image {
    vbox: VoxelBox,
    data: ImageData,
}

impl core::ops::IndexMut<Index1> for Image {
    fn index_mut(&mut self, i: Index1) -> &mut Self::Output { &mut self.data[i] }
}

impl core::ops::Index<Index1> for Image {
    type Output = Intensity;
    fn index(&self, i: Index1) -> &Self::Output { &self.data[i] }
}

#[allow(nonstandard_style)]
impl Image {

    pub fn mlem<'a>(vbox: VoxelBox, measured_lors: &Vec<LOR>/*, noise: &Noise*/) {
        println!("Entered mlem");
        // TODO: sensitivity matrix, all ones for now
        let S = Self::ones(vbox.clone()).data;
        // TODO: noise
        let noise = Noise;
        // TODO: decide how long to iterate
        let however_many = 3;

        // Start off with a uniform image
        let mut image = Self::ones(vbox);

        for n in 0..however_many {
            let BP = image.backproject(measured_lors, &noise);
            for (i, _) in BP.iter().enumerate() { // TODO see if it's any better with iterators
                if S[i] > 0.0 { image[i] *= BP[i] / S[i] }
                else          { image[0]  = 0.0          }
            }
        }
        println!("iteration complete");
    }

    fn backproject<'a>(&'a self, measured_lors: &Vec<LOR>, noise: &Noise) -> ImageData {
        println!("Entered backproject with vbox {:?} and data size {}", self.vbox, self.data.len());
        // TODO: tof sigma
        let sigma = None;

        // Accumulator for all backprojection contributions in this iteration
        let mut BP = self.zeros_buffer();

        // For each measured LOR ...
        for (i, LOR_i) in measured_lors.iter().enumerate() {

            // Weights of all voxels contributing to this LOR
            let A_ijs: Vec<Index1Weight> = LOR_i.active_voxels(&self.vbox, None, sigma).collect();

            // Projection of current image into this LOR
            let P_i = self.project(A_ijs.iter().copied(), noise, i);

            // This LOR's contribution to the backprojection
            for (j, A_ij) in A_ijs {
                BP[j] += A_ij / P_i;
            }
        }

        // Return the total backprojection
        BP
    }

    fn project(&self, A_ijs: impl Iterator<Item = Index1Weight>, b: &Noise, i: usize) -> Weight {
        let lambda = self;
        A_ijs.map(move |(j, A_ij)|  A_ij * lambda[j] + b[i])
            .sum()
    }

    // A new empty data store with matching size
    fn zeros_buffer(&self) -> ImageData { vec![0.0; self.vbox.n_voxels()] }
    fn ones(vbox: VoxelBox) -> Self {
        println!("ones vbox: {:?}", vbox);
        Self { data: vec![1.0; vbox.n_voxels()], vbox}
    }

    fn write_to_persistent_storage(&self) { todo!("Write image to disk") }

    fn iter(&self) -> std::slice::Iter<Intensity> { self.data.iter() }

}

//--------------------------------------------------------------------------------
struct Noise; // TODO

impl core::ops::Index<Index1> for Noise {
    type Output = Intensity;
    fn index(&self, _: Index1) -> &Self::Output {
        // TODO: no noise for now
        &0.0
    }
}


//--------------------------------------------------------------------------------
#[derive(Clone, Debug)]
pub struct VoxelBox {
    pub half_width: Vector,
    pub n: BoxDim,
    pub voxel_size: Vector,
}

impl VoxelBox {

    pub fn new((dx, dy, dz): (Length, Length, Length), (nx, ny, nz): (usize, usize, usize)) -> Self {
        let half_width = Vector::new(dx, dy, dz);
        let n = BoxDim::new(nx, ny, nz);
        let voxel_size =  Self::voxel_size(n, half_width);
            Self { half_width, n, voxel_size, }
    }

    fn voxel_size(n: BoxDim, half_width: Vector) -> Vector {
        // TODO: generalize conversion of VecOf<int> -> VecOf<float>
        let nl: Vector = Vector::new(n.x as Length, n.y as Length, n.z as Length);
        (half_width * 2.0).component_div(&nl)
    }

    pub fn voxel_centre(&self, i: Index3) -> Point {
        //i.map(|n| n as f64 + 0.5).component_mul(&self.voxel_size).into()
        let s = self.voxel_size;
        Point::new((i[0] as f64 + 0.5) * s.x,
                   (i[1] as f64 + 0.5) * s.y,
                   (i[2] as f64 + 0.5) * s.z,)
    }

    pub fn index3_to_1(&self, i: Index3) -> Index1 {
        let n = self.n;
        i[0] + n[0] * i[1] + (n[0] * n[1]) * i[2]
    }

    pub fn index1_to_3(&self, i: Index1) -> Index3 {
        let n = self.n;
        let z = i / (n.x * n.y);
        let r = i % (n.x * n.y);
        let y = r / n.x;
        let x = r % n.x;
        [x,y,z]
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

    fn n_voxels(&self) -> usize {
        let n = &self.n;
        n.x * n.y * n.z
    }
}

#[cfg(test)]
mod test_vbox {
    use super::*;
    use rstest::rstest;

    type I3 = (usize, usize, usize);

    // -------------------- Some hand-picked examples ------------------------------
    #[rstest(/**/    size   , index3 , index1,
             // 1-d examples
             case(( 1, 1, 1), (0,0,0),   0),
             case(( 9, 1, 1), (3,0,0),   3),
             case(( 1, 8, 1), (0,4,0),   4),
             case(( 1, 1, 7), (0,0,5),   5),
             // Counting in binary: note digit reversal
             case(( 2, 2, 2), (0,0,0),   0),
             case(( 2, 2, 2), (1,0,0),   1),
             case(( 2, 2, 2), (0,1,0),   2),
             case(( 2, 2, 2), (1,1,0),   3),
             case(( 2, 2, 2), (0,0,1),   4),
             case(( 2, 2, 2), (1,0,1),   5),
             case(( 2, 2, 2), (0,1,1),   6),
             case(( 2, 2, 2), (1,1,1),   7),
             // Relation to decimal: note reversal
             case((10,10,10), (1,2,3), 321),
             case((10,10,10), (7,9,6), 697),
    )]
    fn hand_picked(size: I3, index3: I3, index1: usize) {
        let irrelevant = (10.0, 10.0, 10.0);
        let vbox = VoxelBox::new(irrelevant, (size.0, size.1, size.2));
        let index3 = [index3.0, index3.1, index3.2];
        assert_eq!(vbox.index3_to_1(index3), index1);
        assert_eq!(vbox.index1_to_3(index1), index3);
    }

    // -------------------- Exhaustive roundtrip testing ------------------------------
    use proptest::prelude::*;

    // A strategy that picks 3-d index limits, and a 1-d index guaranteed to lie
    // within those bounds.
    fn size_and_in_range_index() -> impl Strategy<Value = (I3, usize)> {
        (1..200_usize, 1..200_usize, 1..200_usize)
            .prop_flat_map(|i| (Just(i), 1..(i.0 * i.1 * i.2)))
    }

    proptest! {
        #[test]
        fn index_roundtrip((size, index) in size_and_in_range_index()) {
            let irrelevant = (10.0, 10.0, 10.0);
            let vbox = VoxelBox::new(irrelevant, size);
            let there = vbox.index1_to_3(index);
            let back  = vbox.index3_to_1(there);
            assert_eq!(back, index)
        }

    }
}
//--------------------------------------------------------------------------------

type N = f64;
fn make_gauss(sigma: N) -> impl Fn(N) -> N {
    let root_two_pi = std::f64::consts::PI.sqrt();
    let a = 1.0 / (sigma * root_two_pi);
    move |x| {
        let y = x / sigma;
        let z = y * y;
        a * (-0.5 * z).exp()
    }
}

/// Returns a closure which keeps track of progress along LOR and adjusts voxel
/// weight according to Gaussian TOF distribution. To be mapped over the stream
/// produced by WeightsAlongLOR, as shown in this example:
///
/// ```
/// # use petalo::weights::{VoxelBox, Point, Length, LOR, WeightsAlongLOR, gaussian};
/// // A highly symmetric system for easy testing of features
/// let vbox = VoxelBox::new((30.0, 30.0, 30.0), (5,5,5));
/// let p1 = Point::new(-100.0, 0.0, 0.0);
/// let p2 = Point::new( 100.0, 0.0, 0.0);
/// let (t1, t2) = (12.3, 12.3);
/// let lor = LOR::new(t1, t2, p1, p2);
/// let sigma = 10.0;
///
/// // Generate geometry-dependent voxel weights
/// let active_voxels = WeightsAlongLOR::new(p1, p2, &vbox)
///     // Adjust weights with gaussian TOF factor
///     .map(gaussian(sigma, &lor, &vbox))
///     // Make index more human-friendly (tuple rather than vector)
///     .map(|(i,w)| ((i[0], i[1], i[2]), w))
///     // Store weights in hash map, keyed on voxel index, for easy retrieval
///     .collect::<std::collections::HashMap<(usize, usize, usize), Length>>();
///
/// // Gaussian centred on origin is symmetric about origin
/// assert_eq!(active_voxels.get(&(0,2,2)) , active_voxels.get(&(4,2,2)));
/// assert_eq!(active_voxels.get(&(1,2,2)) , active_voxels.get(&(3,2,2)));
///
/// // Highest probability in the centre
/// assert!   (active_voxels.get(&(0,2,2)) < active_voxels.get(&(1,2,2)));
/// assert!   (active_voxels.get(&(1,2,2)) < active_voxels.get(&(2,2,2)));
/// ```
// TODO: This test is symmetric in t1,t2; need an asymmetric test.
pub fn gaussian(sigma: Length, lor: &LOR, vbox: &VoxelBox) -> Box<dyn FnMut (Index3Weight) -> Index3Weight> {
    let t1_minus_t2 = lor.t1 - lor.t2;
    match vbox.entry(&lor.p1, &lor.p2).map(|ep| (ep-lor.p1).norm()) {
        // If LOR misses the voxel box, we should never receive any voxels
        // weights to adjust.
        None => Box::new(|_| panic!("Cannot adjust for TOF on LOR that misses image volume.")),
        // If LOR does hit some voxels, find the first hit's distance from p1 ...
        Some(p1_to_entry) => {
            // ... and the TOF peak's distance from p1
            let p1_to_peak = 0.5 * ((lor.p1 - lor.p2).norm() + c * (t1_minus_t2));
            // Will keep track of how far we are from the TOF peak
            let mut distance_to_peak = p1_to_peak - p1_to_entry;
            // Specialize gaussian on given sigma
            let gauss = make_gauss(sigma);
            Box::new(move |(index, weight)| {
                // Use the LOR's half-way point in the voxel, for the TOF factor
                let midpoint = distance_to_peak - weight / 2.0;
                // Update distance for next voxel
                distance_to_peak -= weight;
                // Return the weight suppressed by the gaussian TOF factor
                (index, weight * gauss(midpoint))
            })
        },
    }
}


#[derive(Clone, Copy)]
pub struct LOR {
    pub p1: Point,
    pub p2: Point,
    pub t1: Time,
    pub t2: Time,
}

impl LOR {
    pub fn new(t1: Time, t2: Time, p1: Point, p2: Point) -> Self { Self { t1, t2, p1, p2 } }

    pub fn active_voxels<'a>(&self, vbox: &'a VoxelBox, threshold: Option<Length>, tof: Option<Length>) -> impl Iterator<Item = Index1Weight> + 'a {

        //
        let tof_adjustment = match tof {
            Some(sigma) => gaussian(sigma, &self, vbox),
            None        => Box::new(|x| x),
        };

        // Should we ignore voxels with very low contribution?
        let threshold: Box<dyn FnMut(&Index3Weight) -> bool> = match threshold {
            Some(thresh) => Box::new(move |(_, w)| w > &thresh),
            None         => Box::new(     |  _   |  true      ),
        };

        WeightsAlongLOR::new(self.p1, self.p2, vbox)
            .map(tof_adjustment)
            .filter(threshold)
            .map(move |(i, w)| (vbox.index3_to_1(i), w))
    }

}
