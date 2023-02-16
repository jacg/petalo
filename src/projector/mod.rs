pub use siddon::Siddon;

pub mod siddon;

use crate::projector::siddon::{Fs, project_one_lor};
use units::ratio_;

pub fn project_one_lor_mlem<'i, S: SystemMatrix>(fold_state: Fs<'i,S>, lor: &LOR) -> Fs<'i,S> {
    project_one_lor(fold_state, lor, |projection, lor| ratio_(projection * lor.additive_correction))
}

pub fn project_one_lor_sens<'i, S: SystemMatrix>(fold_state: Fs<'i,S>, lor: &LOR) -> Fs<'i,S> {
    project_one_lor(fold_state, lor, |projection, _lor| (-projection).exp())
}


/// Common core of forward and backward propagation.
/// Used by `one_iteration` (MLEM) and `sensitivity_image`
pub fn projector_core<'l, 'i, S, F>(
    projector_data : S,
    image          : &'i Image,
    lors           : &'l [LOR],
    job_size       : usize,
    project_one_lor: F,
) -> ImageData
where
    S: SystemMatrix + Copy + Send + Sync,
    F: Fn(FoldState<'i, S>, &'i LOR) -> FoldState<'i, S> + Sync + Send,
    'l: 'i,
{
    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let initial_thread_state = || {
        let backprojection = Image::zeros_buffer(image.fov);
        let system_matrix_row = S::buffers(image.fov);
        FoldState { backprojection, system_matrix_row, image, projector_data }
    };

    // -------- Project all LORs forwards and backwards ---------------------
    let fold_result = lors
        .par_iter()
        // Rayon is too eager in spawning small jobs, each of which requires the
        // construction and subsequent combination of expensive accumulators
        // (whole `Image`s). So here we try to limit it to one job per thread.
        .with_min_len(job_size)
        .with_max_len(job_size)
        .fold(initial_thread_state, project_one_lor);

    // -------- extract relevant information (backprojection) ---------------
    fold_result
        // Keep only the backprojection (ignore weights and indices)
        .map(|tuple| tuple.backprojection)
        // Sum the backprojections calculated on each thread
        .reduce(|| Image::zeros_buffer(image.fov), elementwise_add)
}

pub fn elementwise_add(a: Vec<f32>, b: Vec<f32>) -> Vec<f32> {
    a.iter().zip(b.iter()).map(|(l,r)| l+r).collect()
}

// Data needed by project_one_lor, both as input and output, because of the
// constrains imposed by `fold`
pub struct FoldState<'img, T: ?Sized> {
    backprojection: ImageData,
    system_matrix_row: SystemMatrixRow,
    image: &'img Image,
    projector_data: T,
}


use rayon::prelude::*;

use crate::{
    LOR,
    image::{ImageData, Image},
    system_matrix::{SystemMatrixRow, SystemMatrix},
};
