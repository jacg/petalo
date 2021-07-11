use crate::types::{Length, TWOPI, ps_to_mm};

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
    sigma.map(|sigma| make_gauss(ps_to_mm(sigma), cutoff))
}
