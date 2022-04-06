use crate::{Length, PerLength, Ratio, Time, TWOPI, C};

use geometry::uom::uomcrate as guomc;
use guomc::ConstZero; // num_traits::Zero;
use geometry::uom::mm;

// How would you make this generic over Length -> T ?
fn make_gauss(sigma: Length, cutoff: Option<Ratio>) -> impl Fn(Length) -> PerLength {
    let two_pi: Ratio = TWOPI;
    let root_two_pi: Ratio = two_pi.sqrt();
    let peak_height: PerLength = 1.0 / (sigma * root_two_pi);
    let cutoff: Length = cutoff.map_or(mm(std::f32::INFINITY), |width| width * sigma);
    move |dx: Length| -> PerLength {
        if dx.abs() < cutoff {
            let y: Ratio = dx / sigma;
            let z: Ratio = y * y;
            peak_height * (-0.5 * z).get::<guomc::si::ratio::ratio>().exp()
        } else {
            PerLength::ZERO
        }
    }
}

pub fn make_gauss_option(sigma: Option<Time>, cutoff: Option<Ratio>) -> Option<impl Fn(Length) -> PerLength> {
    sigma.map(|sigma| make_gauss(sigma * C, cutoff))
}
