use ndarray::azip;

#[cfg(not(feature = "serial"))]
use rayon::prelude::*;

use crate::types::{Length, Time, Ratio};
use crate::weights::{lor_vbox_hit, find_active_voxels, VoxelBox, LOR};
use crate::gauss::make_gauss_option;

type Intensity = Length;

pub type ImageData = Vec<Intensity>;

type Index1 = usize;


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

        // TOF adjustment to apply to the weights
        let tof: Option<_> = make_gauss_option(sigma, cutoff);

        // Closure preparing the state needed by each thread: will be called by
        // `fold` when a thread is launched.
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
            (backprojection, weights, indices)
        };

        // Parallel fold takes a function which will return ID value;
        // serial fold takes the ID value itself.
        #[cfg (feature = "serial")]
        // In the serial case, call the function to get one ID value
        let initial_thread_state =  initial_thread_state();

        // Choose between serial parallel iteration
        #[cfg    (feature = "serial") ] let iter = measured_lors.    iter();
        #[cfg(not(feature = "serial"))] let iter = measured_lors.par_iter();

        // For each measured LOR ...
        let final_thread_state = iter.fold(
            // Empty accumulator (backprojection) and temporary workspace (weights, items)
            initial_thread_state,
            // Process one LOR, storing contribution in `backprojection`
            |(mut backprojection, mut weights, mut indices), lor| {

                // Analyse point where LOR hits voxel box
                match lor_vbox_hit(lor, self.vbox) {

                    // LOR missed voxel box: nothing to be done
                    None => return (backprojection, weights, indices),

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
                // Return the final state collected on this thread
                (backprojection, weights, indices)
            }
        );

        // In the serial case, there is a single result to unwrap ...
        #[cfg (feature = "serial")]
        let backprojection = final_thread_state.0; // Keep only backprojection

        // ... in the parallel case, the results from each thread must be
        // combined
        #[cfg(not(feature = "serial"))]
        let backprojection = {
            final_thread_state
            // Keep only the backprojection (ignore weights and indices)
            .map(|tuple| tuple.0)
            // Sum the backprojections calculated on each thread
            .reduce(|   | self.zeros_buffer(),
                    |l,r| l.iter().zip(r.iter()).map(|(l,r)| l+r).collect())
        };

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
