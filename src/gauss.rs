use crate::types::{Length, Ratio, C, TWOPI};

pub const DISTANCE_AS_TOF_DELTA: Length = 2.0 / C;
pub const TOF_DT_AS_DISTANCE: Ratio = C / 2.0;

fn make_gauss(sigma: Length, cutoff: Option<Length>) -> impl Fn(Length) -> Length {
    let root_two_pi = TWOPI.sqrt() as Length;
    let peak_height = 1.0 / (sigma * root_two_pi);
    let cutoff = cutoff.map_or(std::f32::INFINITY as Length, |width| width * sigma);
    move |x| {
        if x.abs() < cutoff {
            let y = x / sigma;
            let z = y * y;
            peak_height * (-0.5 * z).exp()
        } else {
            0.0
        }
    }
}

pub fn make_gauss_option(sigma: Option<Length>, cutoff: Option<Length>) -> Option<impl Fn(Length) -> Length> {
    sigma.map(|sigma| make_gauss(sigma * TOF_DT_AS_DISTANCE, cutoff))
}
