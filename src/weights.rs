use nalgebra   as na;
use ncollide3d as nc;
use nc::query::RayCast;

// TODO: no thought has been given to what should be public or private.

pub type Length = f64;
type Vector = nc::math ::Vector<Length>;
pub type Point  = nc::math ::Point <Length>;
type Cuboid = nc::shape::Cuboid<Length>;
type Ray    = nc::query::Ray   <Length>;
type Isometry = nc::math::Isometry<Length>;

type VecOf<T> = nc::math::Vector<T>;

pub type Index = VecOf<usize>;
type BoxDim = VecOf<usize>;

/// ```
/// assert_eq!(2+3, 5);
/// ```

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
        index: Index,

        // We exploit the symmetries of the system by flipping some axes. Here
        // we record which axes have been flipped and must be adjusted when
        // calculating indices to be yielded.
        flipped: VecOf<bool>,
    }
}

impl WeightsAlongLOR {
    pub fn new(mut p1: Point, mut p2: Point, vbox: VoxelBox) -> Self {

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
        let lor_direction: Vector = (p2-p1).normalize();
        let mut entry_point: Point = match {
            let lor_length: Length = (p2 - p1).norm();
            let iso: Isometry = Isometry::identity();
            let lor: Ray = Ray::new(p1, lor_direction);
            vbox.aabb.toi_with_ray(&iso, &lor, lor_length, true)
                .map(|toi| lor.origin + lor.dir * toi)
        } { // If LOR misses the box, immediately return an iterator which will
            // generate no hits.
            None => return Self::Outside,
            // Otherwise, unwrap the point and continue calculating a more
            // detailed iterator state
            Some(point) => point,
        };

        // Transform coordinates to align box with axes: making the lower
        // boundaries of the box lie on the zero-planes.
        entry_point += vbox.aabb.half_extents;

        // Express entry point in voxel coordinates: floor(position) = index of voxel.
        let mut entry_point: Vector = entry_point.coords.component_div(&vbox.voxel_size());

        // Floating-point subtractions which should give zero, usually miss very
        // slightly: if this error is negative, the next step (which uses floor)
        // will pick the wrong voxel. Work around this problem by assuming that
        // anything very close to zero is exactly zero.
        entry_point.iter_mut().for_each(|x| if x.abs() < 1e-7 { *x = 0.0 });

        // Find N-dimensional index of voxel at entry point.
        let index: Index = entry_point.map(|x| x.floor() as usize);

        // Voxel size in LOR length units: how far must we move along LOR to
        // traverse one voxel, in any dimension.
        let voxel_size: Vector = vbox.voxel_size().component_div(&lor_direction);

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
    type Item = (Index, Length);

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
                let mut true_index = Index::zeros();
                for n in 0..index.len() {
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

    const TWOPI: Length = std::f64::consts::PI as Length * 2.0;

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
        let vbox = VoxelBox::new((size.0, size.1), (n.0, n.1));

        //crate::visualize::lor_weights(p1, p2, vbox.clone());

        // Collect hits
        let hits: Vec<(Index, Length)> =
            WeightsAlongLOR::new(p1, p2, vbox)
            .inspect(|(is, l)| println!("  ({} {})   {}", is.x, is.y, l))
            .collect();

        // Check total length through voxel box
        let total_length: Length = hits.iter()
            .map(|(_index, weight)| weight)
            .sum();
        assert_approx_eq!(total_length, length);

        // Check voxels hit
        let voxels: Vec<(usize, usize)> = hits.into_iter()
            .map(|(index, _weight)| (index.x, index.y))
            .collect();
        assert_eq!(voxels, expected_voxels)
    }

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
            //p1_z     in -200.0..200.0,
            //p2_z     in -200.0..200.0,
            // Voxel box
            dx in  100.0..(150.0 as Length),
            dy in  100.0..(150.0 as Length), //dz in  100.0..180.0,
            nx in  5..50_usize,
            ny in  5..50_usize,              //nz in  5..50,
        ) {
            let p1_theta: Length = p1_angle * TWOPI;
            let p2_theta: Length = p1_theta + (p2_delta * TWOPI);
            let p1 = Point::new(r * p1_theta.cos(), r * p1_theta.sin(), 0.0);
            let p2 = Point::new(r * p2_theta.cos(), r * p2_theta.sin(), 0.0);
            let vbox = VoxelBox::new((dx, dy), (nx, ny));

            // // Values to plug in to visualizer:
            // println!("let p1 = Point::new({}, {});", p1.x, p1.y);
            // println!("let p2 = Point::new({}, {});", p2.x, p2.y);
            // println!("let vbox = VoxelBox::new(({}, {}), ({}, {}));",
            //          vbox.aabb.half_extents.x, vbox.aabb.half_extents.y, vbox.n.x, vbox.n.y);

            let summed: Length = WeightsAlongLOR::new(p1, p2, vbox.clone())
                .inspect(|(is, l)| println!("  ({} {}) {}", is.x, is.y, l))
                .map(|(_index, weight)| weight)
                .sum();

            let entry_point = |p1: Point, p2:Point| {
                let lor_direction: Vector = (p2-p1).normalize();
                let lor_length: Length = (p2 - p1).norm();
                let iso: Isometry = Isometry::identity();
                let lor: Ray = Ray::new(p1, lor_direction);
                vbox.aabb.toi_with_ray(&iso, &lor, lor_length, true)
                    .map(|toi| lor.origin + lor.dir * toi)
            };

            let a = entry_point(p1, p2);
            let b = entry_point(p2, p1);

            let in_one_go = match (a,b) {
                (Some(a), Some(b)) => (a - b).magnitude(),
                _ => 0.0
            };

            assert_approx_eq!(summed, in_one_go);
        }
    }
}
//--------------------------------------------------------------------------------
#[derive(Clone, Debug)]
pub struct VoxelBox {
    pub aabb: Cuboid, // Axis-Aligned Bounding Box
    pub n: BoxDim,
}

impl VoxelBox {

    pub fn new((dx, dy): (Length, Length), (nx, ny): (usize, usize)) -> Self {
        Self {
            aabb: Cuboid::new(Vector::new(dx, dy, 1.0)),
            n: BoxDim::new(nx, ny, 1)
        }
    }

    pub fn voxel_size(&self) -> Vector {
        // TODO: generalize conversion of VecOf<int> -> VecOf<float>
        let n = &self.n;
        let nl: Vector = Vector::new(n.x as Length, n.y as Length, n.z as Length);
        (self.aabb.half_extents * 2.0).component_div(&nl)
    }

}
