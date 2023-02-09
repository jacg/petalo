use ndarray::azip;

use rayon::prelude::*;

use units::todo::{Lengthf32, Intensityf32};

use crate::{
    LOR,
    config::mlem::Tof,
    fov::FOV,
    gauss::make_gauss_option,
    image::{Image, ImageData},
    projector::Projector,
};

pub static mut N_MLEM_THREADS: usize = 1;

pub fn mlem<P: Projector>(
    fov          : FOV,
    measured_lors: &[LOR],
    tof          : Option<Tof>,
    sensitivity  : Option<Image>,
    n_subsets    : usize,
) -> impl Iterator<Item = (Image, Osem)> + '_ {

    // Start off with a uniform image
    let mut image = Image::ones(fov);

    let sensitivity = sensitivity.or_else(|| Some(Image::ones(fov))).unwrap();

    let mut osem = Osem::new(n_subsets);

    // Return an iterator which generates an infinite sequence of images,
    // each one made by performing one MLEM iteration on the previous one
    std::iter::from_fn(move || {
        one_iteration::<P>(&mut image, osem.subset(measured_lors), &sensitivity.data, tof);
        Some((image.clone(), osem.next())) // TODO see if we can sensibly avoid cloning
    })
}

fn one_iteration<P: Projector>(
    image        : &mut Image,
    measured_lors: &[LOR],
    sensitivity  : &[Intensityf32],
    tof          : Option<Tof>
) {
    // -------- Prepare state required by serial/parallel fold --------------

    // TOF adjustment to apply to the weights
    let tof: Option<_> = make_gauss_option(tof);

    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let previous_image_shared = &*image;
    let initial_thread_state = || {
        let backprojection_image_per_thread = Image::zeros_buffer(image.fov);
        let system_matrix_row = P::buffers(image.fov);
        (backprojection_image_per_thread, system_matrix_row, previous_image_shared, &tof)
    };

    // -------- Project all LORs forwards and backwards ---------------------
    let n_mlem_threads = unsafe {
        // SAFETY: modified only once, at the beginning of bin/mlem.rs::main()
        N_MLEM_THREADS
    };
    let job_size = measured_lors.len() / n_mlem_threads;
    let fold_result = measured_lors
        .par_iter()
        // Rayon is too eager in spawning small jobs, each of which requires the
        // construction and subsequent combination of expensive accumulators
        // (whole `Image`s). So here we try to limit it to one job per thread.
        .with_min_len(job_size)
        .with_max_len(job_size)
        .fold(initial_thread_state, P::project_one_lor);

    // -------- extract relevant information (backprojection) ---------------
    let backprojection = fold_result
        // Keep only the backprojection (ignore weights and indices)
        .map(|tuple| tuple.0)
        // Sum the backprojections calculated on each thread
        .reduce(|| Image::zeros_buffer(image.fov), elementwise_add);

    // -------- Correct for attenuation and detector sensitivity ------------
    apply_sensitivity_image(&mut image.data, &backprojection, sensitivity);
}

pub fn elementwise_add(a: Vec<f32>, b: Vec<f32>) -> Vec<f32> {
    a.iter().zip(b.iter()).map(|(l,r)| l+r).collect()
}

fn apply_sensitivity_image(image: &mut ImageData, backprojection: &[Lengthf32], sensitivity: &[Intensityf32]) {
    //  TODO express with Option<matrix> and mul reciprocal
    // Apply Sensitivity matrix
    azip!((voxel in image, &b in backprojection, &s in sensitivity) {
        if s > 0.0 { *voxel *= b * s }
        else       { *voxel  = 0.0   }
    })
}

#[derive(Debug, Clone, Copy)]
pub struct Osem {
    pub n_subsets: usize,
    pub iteration: usize,
    pub subset: usize,
}

impl Osem {
    fn new(n_subsets: usize) -> Self {
        Self { n_subsets, iteration: 1, subset: 1 }
    }

    fn next(&mut self) -> Self {
        self.subset += 1;
        if self.subset > self.n_subsets {
            self.subset = 1;
            self.iteration += 1;
        }
        *self
    }

    fn subset<'s, 'l>(&'s self, lors: &'l [LOR]) -> &'l [LOR] {
        let set_size = lors.len() / self.n_subsets; // TODO remainder LORs ignored
        let lo = (self.subset - 1) * set_size;
        let hi = lo + set_size;
        &lors[lo..hi]
    }
}
