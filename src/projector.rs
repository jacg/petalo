//! Overall structure of forward and backward projections.
//!
//! The main projection driver functions
//!
//! + `project_lors`, which loops over ...
//!
//! + `project_one_lor`
//!
//! are abstracted over different algorithms for calculating system matrix
//! elements, via the `Projector` trait.
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
/// via the `Projector` trait.
///
/// Different varieties of projection are made possible by injecting the
/// `project_one_lor` function, for which two implementations are provided:
///
/// + `project_one_lor_mlem`
///
/// + `project_one_lor_sens`

use std::borrow::Borrow;
pub fn project_lors<'i, S, L, F>(
    lors          : impl IntoParallelIterator<Item = L>,
    projector_data: S::Data,
    image         : &'i Image,
    result_fov    : Option<FOV>, // if different from image
    project_one_lor: F,
) -> ImageData
where
    S: Projector,
    L: Borrow<LOR>,
    F: Fn(Fs<'i, S>, L) -> Fs<'i, S> + Sync + Send,
{
    let lors = lors.into_par_iter();
    let fwd_fov = image.fov;
    let bck_fov = result_fov.unwrap_or(fwd_fov);
    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let initial_thread_state = || {
        let backprojection = Image::zeros_buffer(result_fov.unwrap_or(bck_fov));
        let matrix_row_fwd =                               S::buffers(bck_fov);
        let bck = result_fov.map(|fov|               (fov, S::buffers(fwd_fov)));
        Fs::<S> { backprojection, matrix_row_fwd, bck, image, projector_data }
    };

    lors
        .fold(initial_thread_state, project_one_lor)
        // Keep only the backprojection (ignore weights and indices)
        .map(|state| state.backprojection)
        // Sum the backprojections calculated on each thread
        .reduce(|| Image::zeros_buffer(bck_fov), elementwise_add)
}

// ----- For injection into `project_lors` --------------------------------------------------
/// Adapts `project_lors` for MLEM iterations
pub fn project_one_lor_mlem<'i, S: Projector>(fold_state: Fs<'i,S>, lor: &LOR) -> Fs<'i,S> {
    project_one_lor::<S>(fold_state, lor, |projection, lor| ratio_(projection * lor.additive_correction))
}

/// Adapts `project_lors` for sensitivity image generation
pub fn project_one_lor_sens<S: Projector>(fold_state: Fs<S>, lor: impl Borrow<LOR>) -> Fs<S> {
    project_one_lor::<S>(fold_state, lor.borrow(), |projection, _lor| (-projection).exp())
}
// ---------------------------------------------------------------------------------------------

/// Used by `project_lors` to perform the forward and backward projection of a
/// single LOR.
///
/// Abstracted over different algorithms for calculating system matrix elements,
/// via the `Projector` trait.
///
/// Different varieties of projection are made possible by injecting the
/// `adapt_forward_projection` function.
fn project_one_lor<'img, S: Projector>(
    state: Fs<'img, S>,
    lor: &LOR,
    adapt_forward_projection: impl Fn(f32, &LOR) -> f32,
) -> Fs<'img, S> {
    let Fs::<S> { mut backprojection, mut matrix_row_fwd, mut bck, image, projector_data } = state;
    matrix_row_fwd.clear(); // Throw away previous LOR's values
    S::update_system_matrix_row(&mut matrix_row_fwd, lor, image.fov, &projector_data);

    // If forward- and back-projection geometries differ, calculate backprojection geometry matrix
    if let Some((fov, matrix_row_bck)) = bck.as_mut() {
        matrix_row_bck.clear(); // Throw away previous LOR's values
        S::update_system_matrix_row(matrix_row_bck, lor, *fov, &projector_data);
    };
    {
        // Does backprojection use same geometry matrix as forward projection?
        let matrix_row_bck = bck.as_ref().map_or(&matrix_row_fwd, |x| &x.1);

        // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
        let project_this_lor = 'safe_lor: {
            for (i, _) in matrix_row_bck {
                if i >= backprojection.len() { break 'safe_lor false; }
            }
            // If bck/fwd projections differ, need a second check
            if bck.is_some() {
                for (i, _) in &matrix_row_fwd {
                    if i >= image.data.len() { break 'safe_lor false; }
                }
            }
            // This LOR looks safe: process it
            true
        };

        if project_this_lor {
            // Sum product of relevant voxels' weights and activities
            let projection = forward_project(&matrix_row_fwd, image);

            // ... the sum needs to be adapted for the specific use case: MLEM
            // or sensitivity image generation, are the only ones so far
            let adapted_projection = adapt_forward_projection(projection, lor);

            // Backprojection of LOR onto image
            back_project(&mut backprojection, matrix_row_bck, adapted_projection);
        }
    }
    // Return values needed by next LOR's iteration
    Fs::<S> { backprojection, matrix_row_fwd, bck, image, projector_data }
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
    pub matrix_row_fwd:   SystemMatrixRow,
    pub bck: Option<(FOV, SystemMatrixRow)>,
    pub image: &'img Image,
    pub projector_data: T,
}

/// Helper for `FoldState`: removes the need to repeat `::Data` at points of
/// use.
pub type Fs<'i, S> = FoldState<'i, <S as Projector>::Data>;

// ----- Imports ------------------------------------------------------------------------------------------
use rayon::prelude::*;

use units::{
    todo::Lengthf32,
    ratio_,
};

use crate::{
    LOR, FOV,
    image::{ImageData, Image},
    projectors::{SystemMatrixRow, Projector},
};
