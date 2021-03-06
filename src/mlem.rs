use ndarray::azip;

#[cfg(not(feature = "serial"))]
use rayon::prelude::*;

use crate::types::{Length, Time, Ratio, Index1};
use crate::weights::{lor_vbox_hit, find_active_voxels, VoxelBox, LOR};
use crate::gauss::make_gauss_option;

type Intensity = Length;

pub type ImageData = Vec<Intensity>;


#[derive(Clone)]
pub struct Image {
    vbox: VoxelBox,
    pub data: ImageData,
    size: usize,
}

impl core::ops::IndexMut<Index1> for Image {
    #[inline]
    fn index_mut(&mut self, i: Index1) -> &mut Self::Output { &mut self.data[i] }
}

impl core::ops::Index<Index1> for Image {
    type Output = Intensity;
    #[inline]
    fn index(&self, i: Index1) -> &Self::Output { &self.data[i] }
}

impl Image {

    pub fn mlem<'a>(vbox: VoxelBox,
                    measured_lors: &'a [LOR],
                    sigma        :     Option<Time>,
                    cutoff       :     Option<Ratio>,
                    smatrix      : &'a [Intensity],
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

    fn one_iteration(&mut self, measured_lors: &[LOR], smatrix: &[Intensity], sigma: Option<Time>, cutoff: Option<Ratio>) {

        // -------- Prepare state required by serial/parallel fold --------------

        // TOF adjustment to apply to the weights
        let tof: Option<_> = make_gauss_option(sigma, cutoff);

        // Closure preparing the state needed by `fold`: will be called by
        // `fold` at the start of every thread that is launched.
        let initial_thread_state = || {

            // Accumulator for all backprojection contributions in this iteration
            let backprojection = self.zeros_buffer();

            // Storage space for the weights and indices of the active voxels
            // (allocating new result vectors for each LOR had a noticeable runtime cost)
            let (weights, indices) = {
                let [nx, ny, nz] = self.vbox.n;
                let max_number_of_active_voxels_possible = nx + ny + nz - 2;
                (Vec::with_capacity(max_number_of_active_voxels_possible),
                 Vec::with_capacity(max_number_of_active_voxels_possible))
            };
            (backprojection, weights, indices, &self, &tof)
        };

        // Parallel fold takes a function which will return ID value;
        // serial fold takes the ID value itself.
        #[cfg (feature = "serial")]
        // In the serial case, call the function to get one ID value
        let initial_thread_state =  initial_thread_state();

        // Choose between serial parallel iteration
        #[cfg    (feature = "serial") ] let iter = measured_lors.    iter();
        #[cfg(not(feature = "serial"))] let iter = measured_lors.par_iter();

        // -------- find active voxels in all LORs ------------------------------

        let fold_result = iter.fold(initial_thread_state, project_one_lor);

        // -------- extract relevant information (backprojection) ---------------

        // In the serial case, there is a single result to unwrap ...
        #[cfg (feature = "serial")]
        let backprojection = fold_result.0; // Keep only backprojection

        // ... in the parallel case, the results from each thread must be
        // combined
        #[cfg(not(feature = "serial"))]
        let backprojection = {
            fold_result
            // Keep only the backprojection (ignore weights and indices)
            .map(|tuple| tuple.0)
            // Sum the backprojections calculated on each thread
            .reduce(|   | self.zeros_buffer(),
                    |l,r| l.iter().zip(r.iter()).map(|(l,r)| l+r).collect())
        };

        // ----------------------------------------------------------------------

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


type FoldState<'r, 'i, 'g, G> = (ImageData , Vec<Length>, Vec<Index1> , &'r &'i mut Image, &'g Option<G>);

fn project_one_lor<'r, 'i, 'g, G>(state: FoldState<'r, 'i, 'g, G>, lor: &LOR) -> FoldState<'r, 'i, 'g, G>
where
    G: Fn(Length) -> Length
{
    let (mut backprojection, mut weights, mut indices, image, tof) = state;

    // Analyse point where LOR hits voxel box
    match lor_vbox_hit(lor, image.vbox) {

        // LOR missed voxel box: nothing to be done
        None => return (backprojection, weights, indices, image, tof),

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
            let projection = forward_project(&weights, &indices, image);

            // Backprojection of LOR onto image
            back_project(&mut backprojection, &weights, &indices, projection);
        }
    }
    // Return the final state collected on this thread
    (backprojection, weights, indices, image, tof)
}

#[inline]
fn forward_project(weights: &[Length], indices: &[usize], image: &Image) -> Length {
    let mut projection = 0.0;
    for (w, j) in weights.iter().zip(indices.iter()) {
        projection += w * image[*j]
    }
    projection
}

#[inline]
fn back_project(backprojection: &mut Vec<Length>, weights: &[Length], indices: &[usize], projection: Length) {
    let projection_reciprocal = 1.0 / projection;
    for (w, j) in weights.iter().zip(indices.iter()) {
        backprojection[*j] += w * projection_reciprocal;
    }
}


fn apply_sensitivity_matrix(image: &mut ImageData, backprojection: &[Length], smatrix: &[Intensity]) {
    //  TODO express with Option<matrix> and mul reciprocal
    // Apply Sensitivity matrix
    azip!((voxel in image, &b in backprojection, &s in smatrix) {
        if s > 0.0 { *voxel *= b / s }
        else       { *voxel  = 0.0   }
    })
}

// --------------------------------------------------------------------------------
//                  Conversion between 1d and 3d indices

use std::ops::{Add, Div, Mul , Rem};

pub fn index3_to_1<T>([ix, iy, iz]: [T; 3], [nx, ny, _nz]: [T; 3]) -> T
where
    T: Mul<Output = T> + Add<Output = T>
{
    ix + (iy + iz * ny) * nx
}

#[allow(clippy::many_single_char_names)]
pub fn index1_to_3<T>(i: T, [nx, ny, _nz]: [T; 3]) -> [T; 3]
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
    use crate::types::Index3;

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
