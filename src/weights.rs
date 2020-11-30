use ncollide3d as nc;
use nc::query::RayCast;

use ndarray::azip;
use ndarray::Array3;

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

/// An iterator which yields the N-dimensional indices of voxels which have been
/// traversed by a LOR, along with the distance the LOR covered inside that
/// voxel.
pub enum WeightsAlongLOR {

    // We are at a point outside of the voxel box: no bookkeeping to be done.
    Outside,

    // We are traversing the voxel box:
    Inside {

        // Distance travelled along LOR since entering the voxel box
        here: Length,

        // Position of the next voxel boundary in any dimension, in LOR distance
        // units from point of entry.

        next_boundary: Vector,

        // Voxel box size, as number of voxels in each dimension
        n_voxels: BoxDim,

        // Dimensions of the voxels expressed in LOR distance units. Used to
        // reset components of `next_boundary` when they reach 0.
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
        entry_point.iter_mut().for_each(|x| if x.abs() < EPS { *x = 0.0 });

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
        let next_boundary: Vector = (Vector::repeat(1.0) - vox_done_fraction)
            .component_mul(&voxel_size);

        // Initial iterator state: the point where LOR enters voxel box (in
        // voxel coordinates), along with bookkeeping information.
        Self::Inside { here: 0.0, next_boundary, voxel_size, n_voxels: vbox.n, flipped, index }
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
            Self::Inside { here, flipped, index, n_voxels, next_boundary, voxel_size } => {

                // Remember index of the voxel we are about to cross (flipped
                // back from our algorithm's internal coordinate system, to the
                // client's original coordinate system).
                let mut true_index = [0; 3];
                for n in 0..3 {
                    if flipped[n] { true_index[n] = n_voxels[n] - 1 - index[n]; }
                    else          { true_index[n] =                   index[n]; }
                }

                // Which boundary will be hit next, and where
                let (dimension, boundary_position) = next_boundary.argmin();

                // The weight of this voxel is the distance from where we
                // entered the voxel to the nearest boundary
                let weight = boundary_position - *here;

                // Move to the boundary of this voxel
                *here = boundary_position;

                // Find the next boundary in this dimension
                next_boundary[dimension] += voxel_size[dimension];

                // Change the index according to the boundary we are crossing
                index[dimension] += 1;

                // If we have traversed the whole voxel box
                if index[dimension] >= n_voxels[dimension] {
                    // no more steps need to be taken after this one.
                    *self = Self::Outside
                }

                // Yield the N-dimensional index of the voxel we have just
                // crossed (expressed in the client's coordinate system), along
                // with the distance that the LOR covered in that voxel.
                Some((true_index, weight))
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
        let hits: Vec<Index3Weight> =
            WeightsAlongLOR::new(p1, p2, &vbox)
            .filter(|(_, w)| w > &0.0)
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

            // Values to plug in to visualizer:
            let lor = LOR::new(0.0, 0.0, p1, p2);
            let command = crate::visualize::vislor_command(&vbox, &lor);
            println!("\nTo visualize this case, run:\n{}\n", command);

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

            #[cfg    (feature = "f64") ] assert_approx_eq!(summed, in_one_go);
            #[cfg(not(feature = "f64"))] assert_approx_eq!(summed, in_one_go, 1e-3);

        }
    }
}
//--------------------------------------------------------------------------------
type Intensity = Length;

type ImageData = Array3<Intensity>;

#[derive(Clone)]
pub struct Image {
    vbox: VoxelBox,
    pub data: ImageData,
}

impl core::ops::IndexMut<Index3> for Image {
    fn index_mut(&mut self, i: Index3) -> &mut Self::Output { &mut self.data[i] }
}

impl core::ops::Index<Index3> for Image {
    type Output = Intensity;
    fn index(&self, i: Index3) -> &Self::Output { &self.data[i] }
}

#[allow(nonstandard_style)]
impl Image {

    pub fn mlem<'a>(vbox: VoxelBox,
                    measured_lors: &'a [LOR],
                    sigma        :     Option<Time>,
                    cutoff       :     Option<Ratio>,
                    S            : &'a ImageData,
    ) -> impl Iterator<Item = Image> + 'a {

        // Start off with a uniform image
        let mut image = Self::ones(vbox);

        // Return an iterator which generates an infinite sequence of images,
        // each one made by performing one MLEM iteration on the previous one
        std::iter::from_fn(move || {
            image.one_iteration(measured_lors, S, sigma, cutoff);
            Some(image.clone()) // TODO see if we can sensibly avoid cloning
        })
    }

    fn one_iteration(&mut self, measured_lors: &[LOR], S: &ImageData, sigma: Option<Time>, cutoff: Option<Ratio>) {

        // Accumulator for all backprojection contributions in this iteration
        let mut BP = self.zeros_buffer();

        let tof = sigma.map(|sigma| make_gauss(sigma * TOF_DT_AS_DISTANCE, cutoff));

        let [nx, ny, nz] = self.vbox.n;
        let max_number_of_active_voxels_possible = nx + ny + nz - 2;
        let mut weights      = Vec::with_capacity(max_number_of_active_voxels_possible);
        let mut true_indices = Vec::with_capacity(max_number_of_active_voxels_possible);

        // For each measured LOR ...
        measured_lors
            .iter()
            .for_each(|LOR_i| {

                // Weights of all voxels contributing to this LOR
                //let A_ijs: Vec<Index3Weight> = LOR_i.active_voxels(&self.vbox, cutoff, sigma).collect();

                // ================================================================================
                // fn new

                let n_voxels = self.vbox.n;
                let mut here = 0.0;

                let mut p1 = LOR_i.p1;
                let mut p2 = LOR_i.p2;

                weights.clear();
                true_indices.clear();

                let dimensions = 3;

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
                let mut entry_point: Point = match self.vbox.entry(&p1, &p2) {
                    // If LOR misses the box, immediately return
                    None => return,
                    // Otherwise, unwrap the point and continue calculating a more
                    // detailed iterator state
                    Some(point) => point,
                };

                // ... and how far the entry point is from the TOF peak
                let t1_minus_t2 = LOR_i.t1 - LOR_i.t2;
                let p1_to_peak = 0.5 * ((p1 - p2).norm() + C * (t1_minus_t2));
                let p1_to_entry = (entry_point - p1).norm();
                let tof_peak = p1_to_peak - p1_to_entry;

                // Transform coordinates to align box with axes: making the lower
                // boundaries of the box lie on the zero-planes.
                entry_point += self.vbox.half_width;

                // Express entry point in voxel coordinates: floor(position) = index of voxel.
                let mut entry_point: Vector = entry_point.coords.component_div(&self.vbox.voxel_size);

                // Floating-point subtractions which should give zero, usually miss very
                // slightly: if this error is negative, the next step (which uses floor)
                // will pick the wrong voxel. Work around this problem by assuming that
                // anything very close to zero is exactly zero.
                entry_point.iter_mut().for_each(|x| if x.abs() < EPS { *x = 0.0 });

                // Find N-dimensional index of voxel at entry point.
                let mut index: Index3 = [entry_point.x.floor() as usize,
                                         entry_point.y.floor() as usize,
                                         entry_point.z.floor() as usize];//entry_point.map(|x| x.floor() as usize);

                // Voxel size in LOR length units: how far must we move along LOR to
                // traverse one voxel, in any dimension.
                let lor_direction = (p2-p1).normalize();
                let voxel_size: Vector = self.vbox.voxel_size.component_div(&lor_direction);

                // What fraction of the voxel has already been traversed at the entry
                // point, along any axis.
                let vox_done_fraction: Vector = entry_point - entry_point.map(|x| x.floor());

                // How far we must travel along LOR before hitting next voxel boundary,
                // in any dimension.
                let mut next_boundary: Vector = (Vector::repeat(1.0) - vox_done_fraction)
                    .component_mul(&voxel_size);

                // ================================================================================
                // fn next
                loop {
                    // Remember index of the voxel we are about to cross (flipped
                    // back from our algorithm's internal coordinate system, to the
                    // client's original coordinate system).
                    let mut true_index = [0; 3];
                    for n in 0..3 {
                        if flipped[n] { true_index[n] = n_voxels[n] - 1 - index[n]; }
                        else          { true_index[n] =                   index[n]; }
                    }

                    // Which boundary will be hit next, and where
                    let (dimension, boundary_position) = next_boundary.argmin();

                    // The weight of this voxel is the distance from where we
                    // entered the voxel to the nearest boundary
                    let mut weight = boundary_position - here;
                    if let Some(gauss) = &tof {
                        weight *= gauss(here - tof_peak);
                    }

                    // Move to the boundary of this voxel
                    here = boundary_position;

                    // Find the next boundary in this dimension
                    next_boundary[dimension] += voxel_size[dimension];

                    // Change the index according to the boundary we are crossing
                    index[dimension] += 1;

                    // If we have traversed the whole voxel box
                    if index[dimension] >= n_voxels[dimension] { break; }

                    // Store the N-dimensional index of the voxel we have just
                    // crossed (expressed in the client's coordinate system), along
                    // with the distance that the LOR covered in that voxel.
                    if weight > 0.0 {
                        true_indices.push(true_index);
                        weights     .push(weight);
                    }
                }
                // ================================================================================

                // Forward projection of current image into this LOR
                // let P_i = self.project(A_ijs.iter().copied());
                let mut projection = 0.0;
                for (w, j) in weights.iter().zip(true_indices.iter()) {
                    projection += w * self[*j]
                }
                let projection_reciprocal = 1.0 / projection;

                // This LOR's contribution to the backprojection
                // for (j, A_ij) in A_ijs {
                //     BP[j] += A_ij / P_i;
                // }
                for (w, j) in weights.iter().zip(true_indices.iter()) {
                    BP[*j] += w * projection_reciprocal;
                }


            });

        // Apply Sensitivity matrix
        azip!((voxel in &mut self.data, &b in &BP, &s in S) {
            if s > 0.0 { *voxel *= b / s }
            else       { *voxel  = 0.0   }
        })

    }

    fn project(&self, A_ijs: impl Iterator<Item = Index3Weight>) -> Weight {
        let lambda = self;
        A_ijs.map(move |(j, A_ij)|  A_ij * lambda[j])
            .sum()
    }

    // A new empty data store with matching size
    fn zeros_buffer(&self) -> ImageData {  Array3::zeros( self.vbox.n  ) }
    pub fn ones(vbox: VoxelBox) -> Self {
        Self { data: Array3::ones( vbox.n  ), vbox}
    }

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

    pub fn index3_to_1(&self, i: Index3) -> Index1 {
        let n = self.n;
        i[0] + n[0] * i[1] + (n[0] * n[1]) * i[2]
    }

    #[allow(clippy::many_single_char_names)]
    pub fn index1_to_3(&self, i: Index1) -> Index3 {
        let n = self.n;
        let z = i / (n[0] * n[1]);
        let r = i % (n[0] * n[1]);
        let y = r / n[0];
        let x = r % n[0];
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

pub fn make_gauss(sigma: Length, cutoff: Option<Length>) -> impl Fn(Length) -> Length {
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

// A doctest, just to remind us that we should have some.

// TODO: This test is symmetric in t1,t2; need an asymmetric test.
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
///     .map(gaussian(sigma, &lor, &vbox, None))
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
pub fn gaussian(sigma: Length, lor: &LOR, vbox: &VoxelBox, cutoff: Option<Length>) -> Box<dyn FnMut (Index3Weight) -> Index3Weight> {
    let t1_minus_t2 = lor.t1 - lor.t2;
    match vbox.entry(&lor.p1, &lor.p2).map(|ep| (ep-lor.p1).norm()) {
        // If LOR misses the voxel box, we should never receive any voxels
        // weights to adjust.
        None => Box::new(|_| panic!("Cannot adjust for TOF on LOR that misses image volume.")),
        // If LOR does hit some voxels, find the first hit's distance from p1 ...
        Some(p1_to_entry) => {
            // ... and the TOF peak's distance from p1
            let p1_to_peak = 0.5 * ((lor.p1 - lor.p2).norm() + C * (t1_minus_t2));
            // Will keep track of how far we are from the TOF peak
            let mut distance_to_peak = p1_to_peak - p1_to_entry;
            // Specialize gaussian on given sigma
            let gauss = make_gauss(sigma, cutoff);
            Box::new(move |(index, weight)| {
                // Use the LOR's half-way point in the voxel, for the TOF factor
                let midpoint = distance_to_peak - weight / 2.0;
                // Update distance for next voxel
                distance_to_peak -= weight;
                // Return the weight suppressed by the gaussian TOF factor
                (index, weight * gauss(midpoint * DISTANCE_AS_TOF_DELTA))
            })
        },
    }
}

#[cfg(test)]
mod test_gaussian {

    use super::*;

    #[test]
    fn peak_shift() {

        // Arbitrarily choose to shift the peak by 40.5 mm from the midpoint
        let delta_x = 40.5;

        // That distance reduces/increases the TOF to the nearer/farther LOR
        // endpoints by delta_t = delta_x / c.
        let delta_t = delta_x / C;

        // As one of the TOFs is reduced by this amount, while the other is
        // increased by the same amount, the total difference between the
        // arrival times, differs by 2 times delta_t.
        let t1 = 1234.5678; // arbitrary choice
        let t2 = t1 + 2.0 * delta_t; // it's the difference that matters

        // LOR parallel to x-axis: all voxels will have identical weights,
        // before taking TOF into consideration.
        let p1 = Point::new(-110.0, 0.5, 0.5);
        let p2 = Point::new( 110.0, 0.5, 0.5);

        // Vbox with width 100 mm, one voxel per mm
        let vbox = VoxelBox::new((100.0, 100.0, 100.0), (100, 100, 100));

        // The peak is at x = -40.5 mm, which has x-index 9
        let expected_voxel_x_index = 9;

        let lor = LOR::new(t1, t2, p1, p2);

        // The value of sigma isn't important here. Use a sharp peak to make
        // debugging output easy to understand.
        let sigma = 0.7;

        let actual_voxel_x_index =
        // Generate geometry-dependent (all equal) voxel weights
            WeightsAlongLOR::new(p1, p2, &vbox)
            // Adjust weights with gaussian TOF factor
            .map(gaussian(sigma, &lor, &vbox, None))
            // Show values: will be hidden if test succeeds
            .inspect(|x| println!("{:?}", x))
            // find the index-and-weight of the peak
            .max_by(|&(_, wa), &(_, wb)|
                    if wa > wb { std::cmp::Ordering::Greater}
                    else       { std::cmp::Ordering::Less }).unwrap()
            // Throw away the weight, keep index
            .0
            // extract x-component from 3d-index
            [0];

        // Values to plug in to visualizer (only show when test fails)
        println!("\nTo visualize this case, run:");
        println!("{}\n", crate::visualize::vislor_command(&vbox, &lor));

        assert_eq!(expected_voxel_x_index, actual_voxel_x_index);

    }

    // TODO: test width of peak

    use proptest::prelude::*;
    use assert_approx_eq::assert_approx_eq;

    proptest! {
        #[test]
        fn gaussian_cutoff(
            sigma in 20.0..(400.0 as Length),
            width in  1.0..(  4.0 as Length)
        ) {
            let gauss = make_gauss(sigma, Some(width));
            let inside  = (1.0 - EPS) * width * sigma;
            let outside = (1.0 + EPS) * width * sigma;
            let pos_in  = gauss( inside);
            let pos_out = gauss( outside);
            let neg_in  = gauss(-inside);
            let neg_out = gauss(-outside);
            assert_eq!(pos_out, 0.0);
            assert_eq!(neg_out, 0.0);
            assert_approx_eq!(pos_in, neg_in);
            assert!(pos_in > 0.0);
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct LOR {
    pub p1: Point,
    pub p2: Point,
    pub t1: Time,
    pub t2: Time,
}

impl LOR {
    pub fn new(t1: Time, t2: Time, p1: Point, p2: Point) -> Self { Self { t1, t2, p1, p2 } }

    pub fn active_voxels<'a>(&self, vbox: &'a VoxelBox, cutoff: Option<Length>, tof: Option<Length>) -> impl Iterator<Item = Index3Weight> + 'a {

        //
        let tof_adjustment = match tof {
            Some(sigma) => gaussian(sigma, &self, vbox, cutoff),
            None        => Box::new(|x| x),
        };

        WeightsAlongLOR::new(self.p1, self.p2, vbox)
            .map(tof_adjustment)
            // Gaussian-truncation will have set some of the voxels' weights to
            // zero. They will contribute nothing to the calculation, so keeping them
            // will merely waste CPU cycles: filter them out
            .filter(|(_, w)| w > &0.0)
    }

}
