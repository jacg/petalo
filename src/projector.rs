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
pub fn project_lors<'i, S, L, F>(
    projector_data : S::Data,
    image          : &'i Image,
    lors           : L,
    job_size       : usize,
    project_one_lor: F,
) -> ImageData
where
    S: SystemMatrix,
    L: IntoParallelIterator,
    L::Iter: IndexedParallelIterator,
    F: Fn(Fs<'i, S>, L::Item) -> Fs<'i, S> + Sync + Send,
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
        .fold_chunks(job_size, initial_thread_state, project_one_lor);

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


// ----- WIP ----------------------------------------------------------------------------------------------
use units::{
    Length,
    mm, ns, ratio
};

use crate::utils::group_digits;

fn WIP_make_points<S>(
    detector_length  : Length,
    detector_diameter: Length,
) -> Vec<crate::Point>
where
    S: SystemMatrix,
{
    dbg!(detector_length);
    // For prototyping purposes, hard-wire the scintillator element size
    let dz = mm(3.0);
    let da = mm(3.0);
    let dr = mm(30.0);
    // Points at the centres of all elements
    let points = crate::discrete::Discretize::new(detector_diameter, dr, dz, da)
        .all_element_centres(detector_length)
        .collect::<Vec<_>>();
    dbg!(group_digits(points.len()));
    points
}

fn WIP_make_lors_par_iter<S>(points: &[crate::Point], fov: crate::FOV) -> impl ParallelIterator<Item = LOR> + '_
where
    S: SystemMatrix,
{
    (0..points.len())
        .par_bridge()
        .flat_map_iter(|i| (i..points.len()).zip(std::iter::repeat(i)))
        .map   (move | (i,j)| (points[i], points[j]))
        .filter(move |&(p,q)| fov.entry(p,q).is_some())
        .map   (move | (p,q)| LOR::new(ns(0.0), ns(0.0), p, q, ratio(1.0)))
}

fn WIP_project_lors<S: SystemMatrix>(
    lors          : impl ParallelIterator<Item = LOR>,
    projector_data: S::Data,
    image         : &Image,
) -> Image {
    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let initial_thread_state = || {
        let backprojection = Image::zeros_buffer(image.fov);
        let system_matrix_row = S::buffers(image.fov);
        Fs::<S> { backprojection, system_matrix_row, image, projector_data }
    };

    let image_data = lors
        .fold(initial_thread_state, |s, l| project_one_lor_sens::<S>(s, &l))
        .map(|state| state.backprojection)
        .reduce(|| Image::zeros_buffer(image.fov), elementwise_add);

    Image::new(image.fov, image_data)
}

pub fn WIP<S>(
    detector_length  : Length,
    detector_diameter: Length,
    projector_data   : S::Data,
    image            : &Image,
) -> Image
where
    S: SystemMatrix,
{
    let points = WIP_make_points       ::<S>(detector_length, detector_diameter);
    let lors   = WIP_make_lors_par_iter::<S>(&points, image.fov);
    let image  = WIP_project_lors      ::<S>(lors, projector_data, image);
    image
}
