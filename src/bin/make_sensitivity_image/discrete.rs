/// Create sensitivity image by backprojecting all possible LORs between pairs
/// of detector elements.
pub fn sensitivity_image<S>(
    detector_length  : Length,
    detector_diameter: Length,
    projector_data   : S::Data,
    attenuation      : &Image,
) -> Image
where
    S: SystemMatrix,
{
    let points     = make_points::<S>(detector_length, detector_diameter);
    let lors       = make_lors  ::<S>(&points, attenuation.fov);
    let mut image_data = project_lors::<S,_,_>(lors, projector_data, attenuation, project_one_lor_sens::<S>);
    let n_lors = { // TODO this is incorrect: it doesn't take the FOV filter into account
        let n = points.len();
        n * (n-1) / 2
    };
    normalize(&mut image_data, n_lors);
    Image::new(attenuation.fov, image_data)
}

fn make_points<S>(
    detector_length  : Length,
    detector_diameter: Length,
) -> Vec<Point>
where
    S: SystemMatrix,
{
    dbg!(detector_length);
    // For prototyping purposes, hard-wire the scintillator element size
    let dz = mm(3.0);
    let da = mm(3.0);
    let dr = mm(30.0);
    // Points at the centres of all elements
    let points = petalo::discrete::Discretize::new(detector_diameter, dr, dz, da)
        .all_element_centres(detector_length)
        .collect::<Vec<_>>();
    dbg!(petalo::utils::group_digits(points.len()));
    points
}

pub (crate) fn make_lors<S>(points: &[Point], fov: crate::FOV) -> impl ParallelIterator<Item = LOR> + '_
where
    S: SystemMatrix,
{
    // let origin = petalo::Point::new(mm(0.0), mm(0.0), mm(0.0));
    (0..points.len())
        .par_bridge()
        .flat_map_iter(|i| (i..points.len()).zip(std::iter::repeat(i)))
        .map   (move | (i,j)| (points[i], points[j]))
        // Rough approximation to 'passes through FOV'
        //.filter(|&q| origin.distance_to_line(*p, *q) < fov.half_width.z)
        .filter(move |&(p,q)| fov.entry(p,q).is_some())
        .map   (move | (p,q)| LOR::new(ns(0.0), ns(0.0), p, q, ratio(1.0)))
}


// ----- Imports ------------------------------------------------------------------------------------------

use petalo::{
    LOR, Point,
    image::Image,
    projector::{project_lors, project_one_lor_sens},
    system_matrix::SystemMatrix,
};

use units::{
    Length,
    mm, ns, ratio,
};

use super::normalize;
use rayon::prelude::*;
