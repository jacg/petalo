use crate::types::{Length, UomPerLength, UomRatio, UomTime, UOM_TWOPI, UOM_C};

use geometry::uom::uomcrate as guomc;
use guomc::ConstZero; // num_traits::Zero;
use geometry::uom::mm;

// How would you make this generic over Length -> T ?
fn make_gauss(sigma: Length, cutoff: Option<UomRatio>) -> impl Fn(Length) -> UomPerLength {
    let two_pi: UomRatio = UOM_TWOPI;
    let root_two_pi: UomRatio = two_pi.sqrt();
    let peak_height: UomPerLength = 1.0 / (sigma * root_two_pi);
    let cutoff: Length = cutoff.map_or(mm(std::f32::INFINITY), |width| width * sigma);
    move |dx: Length| -> UomPerLength {
        if dx.abs() < cutoff {
            let y: UomRatio = dx / sigma;
            let z: UomRatio = y * y;
            peak_height * (-0.5 * z).get::<guomc::si::ratio::ratio>().exp()
        } else {
            UomPerLength::ZERO
        }
    }
}

pub fn make_gauss_option(sigma: Option<UomTime>, cutoff: Option<UomRatio>) -> Option<impl Fn(Length) -> UomPerLength> {
    sigma.map(|sigma| make_gauss(sigma * UOM_C, cutoff))
}
