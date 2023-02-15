pub use siddon::Siddon;

pub mod siddon;

/// Abstract interface for forward-backward projection implementations
pub trait Projector {
    fn project_one_lor    <'i>(fold_state: FoldState<'i, Self>, lor: &LOR) -> FoldState<'i, Self>;
    fn sensitivity_one_lor<'i>(fold_state: FoldState<'i, Self>, lor: &LOR) -> FoldState<'i, Self>;

    // Sparse storage of the slice through the system matrix which corresponds
    // to the current LOR. Allocating these anew for each LOR had a noticeable
    // runtime cost, so we create them up-front and reuse them.
    // This should probably have a default implementation
    fn buffers(fov: FOV) -> SystemMatrixRow;
}


/// Common core of forward and backward propagation.
/// Used by `one_iteration` (MLEM) and `sensitivity_image`
pub fn projector_core<'l, 'i, P, F>(
    projector_data : P,
    image          : &'i Image,
    lors           : &'l [LOR],
    job_size       : usize,
    project_one_lor: F,
) -> ImageData
where
    P: Projector + Copy + Send + Sync,
    F: Fn(FoldState<'i, P>, &'i LOR) -> FoldState<'i, P> + Sync + Send,
    'l: 'i,
{
    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let initial_thread_state = || {
        let backprojection = Image::zeros_buffer(image.fov);
        let system_matrix_row = P::buffers(image.fov);
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
//pub type FoldState<'img, T> = (ImageData, SystemMatrixRow, &'img Image, T);
pub type FoldState<'img, T> = FldStt<'img, T>;

pub struct FldStt<'img, T: ?Sized> {
    backprojection: ImageData,
    system_matrix_row: SystemMatrixRow,
    image: &'img Image,
    projector_data: T,
}


use rayon::prelude::*;

use crate::{
    LOR,
    fov::FOV,
    image::{ImageData, Image},
    system_matrix::SystemMatrixRow,
};
