/// Create sensitivity image by backprojecting all possible LORs between pairs
/// of detector elements.
pub fn sensitivity_image<S>(
    detector_length  : Length,
    projector_data   : S::Data,
    attenuation      : &Image,
    discretize       : Discretize,
    backproj_fov     : Option<FOV>,
) -> Image
where
    S: SystemMatrix,
{
    let points     = make_points::<S>(detector_length, discretize);
    let lors       = make_lors  ::<S>(&points, attenuation.fov, discretize);
    let mut image_data = project_lors::<S,_,_>(lors, projector_data, attenuation, backproj_fov, project_one_lor_sens::<S,_>(backproj_fov));
    let n_lors = { // TODO this is incorrect: it doesn't take the FOV filter into account
        let n = points.len();
        n * (n-1) / 2
    };
    normalize(&mut image_data, n_lors);
    Image::new(attenuation.fov, image_data)
}

fn make_points<S>(
    detector_length  : Length,
    discretize       : Discretize,
) -> Vec<Point>
where
    S: SystemMatrix,
{
    dbg!(detector_length);
    // Points at the centres of all elements
    let points = discretize
        .centre_all_elements(detector_length)
        .collect::<Vec<_>>();
    dbg!(petalo::utils::group_digits(points.len()));
    points
}

pub (crate) fn make_lors<S>(points: &[Point], fov: crate::FOV, discretize: Discretize) -> impl ParallelIterator<Item = LOR> + '_
where
    S: SystemMatrix,
{
    let adjust = discretize.make_adjust_fn();
    let smear: Box<dyn Fn(Point) -> Point + Sync + Send> = {
        Box::new(move |p: Point|  {
            let (x,y,z) = adjust((p.x, p.y, p.z));
            Point { x, y, z }
        })
    };

    // let origin = petalo::Point::new(mm(0.0), mm(0.0), mm(0.0));
    (0..points.len())
        .par_bridge()
        .flat_map_iter(|i| (i..points.len()).zip(std::iter::repeat(i)))
        .map   (move | (i,j)| (points[i], points[j]))
        .map   (move | (p,q)| ( smear(p),  smear(q)))
        // Rough approximation to 'passes through FOV'
        //.filter(|&q| origin.distance_to_line(*p, *q) < fov.half_width.z)
        .filter(move |&(p,q)| fov.entry(p,q).is_some())
        .map   (move | (p,q)| LOR::new(ns(0.0), ns(0.0), p, q, ratio(1.0)))
}


// ----- Imports ------------------------------------------------------------------------------------------

use petalo::{
    FOV, LOR, Point,
    image::Image,
    projector::{project_lors, project_one_lor_sens},
    system_matrix::SystemMatrix,
    discrete::Discretize,
};

use units::{
    Length,
    ns, ratio,
};

use super::normalize;
use rayon::prelude::*;
