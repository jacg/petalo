//! Overall structure of forward and backward projections.
//!
//! The main projection driver functions
//!
//! + `project_lors`, which loops over ...
//!
//! + `project_one_lor`
//!
//! are abstracted over different algorithms for calculating system matrix
//! elements, via the `SystemMatrix` trait.
//!
//! Projections adapted for different use cases can be created by passing
//! adapter functions into `project_lors` and `project_one_lor`. Two varieties
//! are implemented:
//!
//! + MLEM
//!
//! + Sensitivity image construction

/// Performs forward and backward projections over a collection of LORs.
///
/// Abstracted over different algorithms for calculating system matrix elements,
/// via the `SystemMatrix` trait.
///
/// Different varieties of projection are made possible by injecting the
/// `project_one_lor` function, for which two implementations are provided:
///
/// + `project_one_lor_mlem`
///
/// + `project_one_lor_sens`
pub fn project_lors<'l, 'i, S, F>(
    projector_data : S::Data,
    image          : &'i Image,
    lors           : &'l [LOR],
    job_size       : usize,
    project_one_lor: F,
) -> ImageData
where
    S: SystemMatrix,
    F: Fn(Fs<'i, S>, &'i LOR) -> Fs<'i, S> + Sync + Send,
    'l: 'i,
{
    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let initial_thread_state = || {
        let backprojection = Image::zeros_buffer(image.fov);
        let system_matrix_row = S::buffers(image.fov);
        Fs::<S> { backprojection, system_matrix_row, image, projector_data }
    };

    // -------- Project all LORs forwards and backwards ---------------------
    let fold_result = lors
        .into_par_iter()
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

// ----- For injection into `project_lors` --------------------------------------------------
/// Adapts `project_lors` for MLEM iterations
pub fn project_one_lor_mlem<'i, S: SystemMatrix>(fold_state: Fs<'i,S>, lor: &LOR) -> Fs<'i,S> {
    project_one_lor::<S>(fold_state, lor, |projection, lor| ratio_(projection * lor.additive_correction))
}

/// Adapts `project_lors` for sensitivity image generation
pub fn project_one_lor_sens<'i, S: SystemMatrix>(fold_state: Fs<'i,S>, lor: &LOR) -> Fs<'i,S> {
    project_one_lor::<S>(fold_state, lor, |projection, _lor| (-projection).exp())
}
// ---------------------------------------------------------------------------------------------

/// Used by `project_lors` to perform the forward and backward projection of a
/// single LOR.
///
/// Abstracted over different algorithms for calculating system matrix elements,
/// via the `SystemMatrix` trait.
///
/// Different varieties of projection are made possible by injecting the
/// `adapt_forward_projection` function.
fn project_one_lor<'img, S: SystemMatrix>(
    state: Fs<'img, S>,
    lor: &LOR,
    adapt_forward_projection: impl Fn(f32, &LOR) -> f32,
) -> Fs<'img, S> {
    let Fs::<S> { mut backprojection, mut system_matrix_row, image, projector_data } = state;
    // Throw away previous LOR's values
    system_matrix_row.clear();

    S::update_system_matrix_row(&mut system_matrix_row, lor, image.fov, &projector_data);

    let project_this_lor = 'safe_lor: {
            // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
            for (i, _) in &system_matrix_row {
                if i >= backprojection.len() { break 'safe_lor false; }
            }
            // This LOR looks safe: process it
            true
    };

    if project_this_lor {
        // Sum product of relevant voxels' weights and activities
        let projection = forward_project(&system_matrix_row, image);

        // ... the sum needs to be adapted for the use specific use case:
        // MLEM or sensitivity image generation, are the only ones so far
        let adapted_projection = adapt_forward_projection(projection, lor);

        // Backprojection of LOR onto image
        back_project(&mut backprojection, &system_matrix_row, adapted_projection);
    }
    // Return values needed by next LOR's iteration
    Fs::<S> { backprojection, system_matrix_row, image, projector_data }
}

#[inline]
fn forward_project(system_matrix_row: &SystemMatrixRow, image: &Image) -> Lengthf32 {
    let mut projection = 0.0;
    for (j, w) in system_matrix_row {
        projection += w * image[j]
    }
    projection
}

#[inline]
fn back_project(backprojection: &mut [Lengthf32], system_matrix_row: &SystemMatrixRow, projection: Lengthf32) {
    let projection_reciprocal = 1.0 / projection;
    for (j, w) in system_matrix_row {
        backprojection[j] += w * projection_reciprocal;
    }
}

fn elementwise_add(a: Vec<f32>, b: Vec<f32>) -> Vec<f32> {
    a.iter().zip(b.iter()).map(|(l,r)| l+r).collect()
}

/// Data needed to be passed efficiently between the projection of one LOR and
/// the next. Needs to work in conjunction with `rayon`'s `fold`s.
pub struct FoldState<'img, T: ?Sized> {
    pub backprojection: ImageData,
    pub system_matrix_row: SystemMatrixRow,
    pub image: &'img Image,
    pub projector_data: T,
}

/// Helper for `FoldState`: removes the need to repeat `::Data` at points of
/// use.
pub type Fs<'i, S> = FoldState<'i, <S as SystemMatrix>::Data>;

// ----- Imports ------------------------------------------------------------------------------------------
use rayon::prelude::*;

use units::{
    todo::Lengthf32,
    ratio_,
};

use crate::{
    LOR,
    image::{ImageData, Image},
    system_matrix::{SystemMatrixRow, SystemMatrix},
};
