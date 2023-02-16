use ndarray::azip;

use units::todo::{Lengthf32, Intensityf32};

use crate::{
    LOR,
    fov::FOV,
    image::{Image, ImageData},
    projector::{projector_core, project_one_lor_mlem},
};

pub static mut N_MLEM_THREADS: usize = 1;

pub fn mlem<'a, P: SystemMatrix + Copy + Send + Sync + 'a>(
    projector    : P,
    fov          : FOV,
    measured_lors: &'a [LOR],
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
        one_iteration::<P>(projector, &mut image, osem.subset(measured_lors), &sensitivity.data);
        Some((image.clone(), osem.next())) // TODO see if we can sensibly avoid cloning
    })
}

use crate::system_matrix::SystemMatrix;
fn one_iteration<S: SystemMatrix + Copy + Send + Sync>(
    projector    : S,
    image        : &mut Image,
    measured_lors: &[LOR],
    sensitivity  : &[Intensityf32],
) {
    let n_mlem_threads = unsafe {
        // SAFETY: modified only once, at the beginning of bin/mlem.rs::main()
        N_MLEM_THREADS
    };
    let job_size = measured_lors.len() / n_mlem_threads;
    let backprojection = projector_core(projector, &*image, measured_lors, job_size, project_one_lor_mlem);

    // -------- Correct for attenuation and detector sensitivity ------------
    apply_sensitivity_image(&mut image.data, &backprojection, sensitivity);
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
