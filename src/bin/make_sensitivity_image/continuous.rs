/// Create sensitivity image by backprojecting `n_lors` randomly-generated LORs
/// through `attenuation` image.
pub fn sensitivity_image<'l, S: SystemMatrix>(
    detector_length  : Length,
    detector_diameter: Length,
    parameters       : S::Data,
    attenuation      : &Image,
    n_lors           : usize,
) -> Image {
    let lors = find_lors(n_lors, attenuation.fov, detector_length, detector_diameter).
        into_par_iter();

    let mut backprojection = project_lors::<S,_,_>(lors.clone(), parameters, attenuation, project_one_lor_sens::<S>);
    normalize(&mut backprojection, lors.len());
    Image::new(attenuation.fov, backprojection)
}

/// Return a vector of randomly-generated LORs with endpoints on cylinder, passing through the FOV.
pub fn find_lors(n_lors: usize, fov: FOV, detector_length: Length, detector_diameter: Length) -> Vec<LOR> {
    let (l,r) = (detector_length, detector_diameter / 2.0);
    let one_useful_random_lor = move |_lor_number| {
        loop {
            let p1 = random_point_on_cylinder(l, r);
            let p2 = random_point_on_cylinder(l, r);
            if fov.entry(p1, p2).is_some() {
                return LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0))
            }
        }
    };

    use rayon::prelude::*;
    (0..n_lors)
        .into_par_iter()
        .map(one_useful_random_lor)
        .collect()
}

fn random_point_on_cylinder(l: Length, r: Length) -> petalo::Point {
    use std::f32::consts::TAU;
    use rand::random;
    let z     = l   * (random::<Lengthf32>() - 0.5);
    let theta = TAU *  random::<Lengthf32>();
    let x = r * theta.cos();
    let y = r * theta.sin();
    petalo::Point::new(x, y, z)
}

// ----- Imports ------------------------------------------------------------------------------------------

use petalo::{
    FOV, LOR,
    image::Image,
    projector::{project_lors, project_one_lor_sens},
    system_matrix::SystemMatrix,
};

use units::{
    Length, Time,
    todo::Lengthf32,
    ratio,
    uom::ConstZero,
};

use super::normalize;
use rayon::prelude::*;
