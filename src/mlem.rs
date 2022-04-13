use std::path::Path;
use ndarray::azip;

#[cfg(not(feature = "serial"))]
use rayon::prelude::*;

use crate::{io, Lengthf32, Index1_u, Intensityf32};
use crate::{Length, PerLength, Ratio, Time};
use crate::{fov::{lor_fov_hit, FovHit}, system_matrix::{system_matrix_elements, LOR}};
use crate::fov::FOV;
use crate::gauss::make_gauss_option;
use geometry::uom::ratio_;

use crate::image::{Image, ImageData};

impl Image {

    pub fn mlem<'a>(fov: FOV,
                    measured_lors: &'a [LOR],
                    sigma        :     Option<Time>,
                    cutoff       :     Option<Ratio>,
                    sensitivity  :     Option<Self>,
    ) -> impl Iterator<Item = Image> + '_ {

        // Start off with a uniform image
        let mut image = Self::ones(fov);

        let sensitivity = sensitivity.or_else(|| Some(Self::ones(fov))).unwrap();

        // Return an iterator which generates an infinite sequence of images,
        // each one made by performing one MLEM iteration on the previous one
        std::iter::from_fn(move || {
            image.one_iteration(measured_lors, &sensitivity.data, sigma, cutoff);
            Some(image.clone()) // TODO see if we can sensibly avoid cloning
        })
    }

    pub fn from_raw_file(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        Ok((&crate::io::raw::Image3D::read_from_file(path)?).into())
    }

    pub fn write_to_raw_file(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        io::raw::Image3D::from(self).write_to_file(path)?;
        Ok(())
    }

    // Too much copy-paste code reuse from project_one_lor. This is because the
    // latter (and the functions it uses) was heavily optimized, at the cost of
    // ease of reuse.

    // TODO turn this into a method?
    /// Create sensitivity image by backprojecting LORs. In theory this should
    /// use *all* possible LORs. In practice use a representative sample.
    pub fn sensitivity_image(fov: FOV, density: Self, lors: impl Iterator<Item = crate::system_matrix::LOR>, n_lors: usize, stradivarius: f32) -> Self {
        let a = &fov;
        let b = &density.fov;
        if a.n != b.n || a.half_width != b.half_width {
            panic!("For now, attenuation and output image dimensions must match exactly.")
        }
        // TODO convert from density to attenuation coefficient
        let attenuation = density;
        let (mut image, mut weights, mut indices) = projection_buffers(fov);

        // TOF should not be used as LOR attenuation is independent of decay point
        let notof = make_gauss_option(None, None);

        'lor: for lor in lors {
            // Find active voxels (slice of system matrix) WITHOUT TOF
            // Analyse point where LOR hits FOV
            match lor_fov_hit(&lor, fov) {

                // LOR missed FOV: nothing to be done
                None => continue,

                // Data needed by `system_matrix_elements`
                Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

                    // Throw away previous LOR's values
                    weights.clear();
                    indices.clear();

                    // Find active voxels and their weights
                    system_matrix_elements(
                        &mut indices, &mut weights,
                        next_boundary, voxel_size,
                        index, delta_index, remaining,
                        tof_peak, &notof
                    );

                    // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
                    for i in &indices {
                        if *i >= image.len() { continue 'lor; }
                    }

                    let integral = forward_project(&weights, &indices, &attenuation) / stradivarius;
                    let attenuation_factor = (-integral).exp();
                    //println!("{:<8.2e}  ---  {:<8.2e}", integral, attenuation_factor);
                    // Backprojection of LOR onto image
                    back_project(&mut image, &weights, &indices, attenuation_factor);
                }
            }
        }
        // TODO: Just trying an ugly hack for normalizing the image. Do something sensible instead!
        let size = n_lors as f32;
        for e in image.iter_mut() {
            *e /= size
        }
        Self::new(fov, image)
    }

    fn one_iteration(&mut self, measured_lors: &[LOR], sensitivity: &[Intensityf32], sigma: Option<Time>, cutoff: Option<Ratio>) {

        // -------- Prepare state required by serial/parallel fold --------------

        // TOF adjustment to apply to the weights
        let tof: Option<_> = make_gauss_option(sigma, cutoff);

        // Closure preparing the state needed by `fold`: will be called by
        // `fold` at the start of every thread that is launched.
        let initial_thread_state = || {
            let (backprojection, weights, indices) = projection_buffers(self.fov);
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

        // -------- Project all LORs forwards and backwards ---------------------

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
            .reduce(|   | zeros_buffer(self.fov),
                    |l,r| l.iter().zip(r.iter()).map(|(l,r)| l+r).collect())
        };

        // -------- Correct for attenuation and detector sensitivity ------------

        apply_sensitivity_image(&mut self.data, &backprojection, sensitivity);

    }

    pub fn ones(fov: FOV) -> Self {
        let [x,y,z] = fov.n;
        let size = x * y * z;
        Self { data: vec![1.0; size], fov}
    }

    pub fn new(fov: FOV, data: ImageData) -> Self {
        let [x, y, z] = fov.n;
        if data.len() != x * y * z {
            // TODO change panic to Option or Result
            panic!("Image data does not match dimensions {:?}", fov.n);
        };
        Image { fov, data }
    }

    pub fn empty(fov: FOV) -> Self {
        let [x,y,z] = fov.n;
        Self::new(fov, vec![0.0; x*y*z])
    }

    pub fn inverted(&self) -> Self {
        let mut inverted = self.clone();
        for e in inverted.data.iter_mut() { *e = 1.0 / *e }
        inverted
    }
}

fn projection_buffers(fov: FOV) -> (ImageData, Vec<Lengthf32>, Vec<usize>) {
    // The backprojection (or sensitivity image) being constructed in a
    // given MLEM iteration (or sensitivity image calculation).
    let image = zeros_buffer(fov);
    // Weights and indices are sparse storage of the slice through the
    // system matrix which corresponds to the current LOR. (Allocating these
    // anew for each LOR had a noticeable runtime cost.)
    let [nx, ny, nz] = fov.n;
    let max_number_of_active_voxels_possible = nx + ny + nz - 2;
    let weights = Vec::with_capacity(max_number_of_active_voxels_possible);
    let indices = Vec::with_capacity(max_number_of_active_voxels_possible);
    (image, weights, indices)
}

// A new empty data store with matching size
fn zeros_buffer(fov: FOV) -> ImageData { let [x,y,z] = fov.n; vec![0.0; x*y*z] }


type FoldState<'r, 'i, 'g, G> = (ImageData , Vec<Lengthf32>, Vec<Index1_u> , &'r &'i mut Image, &'g Option<G>);

fn project_one_lor<'r, 'i, 'g, G>(state: FoldState<'r, 'i, 'g, G>, lor: &LOR) -> FoldState<'r, 'i, 'g, G>
where
    G: Fn(Length) -> PerLength
{
    let (mut backprojection, mut weights, mut indices, image, tof) = state;

    // Analyse point where LOR hits FOV
    match lor_fov_hit(lor, image.fov) {

        // LOR missed FOV: nothing to be done
        None => return (backprojection, weights, indices, image, tof),

        // Data needed by `system_matrix_elements`
        Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

            // Throw away previous LOR's values
            weights.clear();
            indices.clear();

            // Find active voxels and their weights
            system_matrix_elements(
                &mut indices, &mut weights,
                next_boundary, voxel_size,
                index, delta_index, remaining,
                tof_peak, tof
            );

            // Forward projection of current image into this LOR
            let projection = forward_project(&weights, &indices, image) * lor.additive_correction;

            // Backprojection of LOR onto image
            back_project(&mut backprojection, &weights, &indices, ratio_(projection));
        }
    }
    // Return updated FoldState
    (backprojection, weights, indices, image, tof)
}

#[inline]
fn forward_project(weights: &[Lengthf32], indices: &[usize], image: &Image) -> Lengthf32 {
    let mut projection = 0.0;
    for (w, &j) in weights.iter().zip(indices.iter()) {
        projection += w * image[j]
    }
    projection
}

#[inline]
fn back_project(backprojection: &mut Vec<Lengthf32>, weights: &[Lengthf32], indices: &[usize], projection: Lengthf32) {
    let projection_reciprocal = 1.0 / projection;
    for (w, &j) in weights.iter().zip(indices.iter()) {
        backprojection[j] += w * projection_reciprocal;
    }
}

fn apply_sensitivity_image(image: &mut ImageData, backprojection: &[Lengthf32], sensitivity: &[Intensityf32]) {
    //  TODO express with Option<matrix> and mul reciprocal
    // Apply Sensitivity matrix
    azip!((voxel in image, &b in backprojection, &s in sensitivity) {
        if s > 0.0 { *voxel *= b * s }
        else       { *voxel  = 0.0   }
    })
}

