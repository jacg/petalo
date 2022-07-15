use crate::{Angle, Length, PerLength, Ratio, TWOPI, C};
use crate::config::mlem::Tof;

use geometry::uom::ConstZero; // num_traits::Zero;
use geometry::units::mm;

// How would you make this generic over Length -> T ?
fn make_gauss(sigma: Length, cutoff: Option<Ratio>) -> impl Fn(Length) -> PerLength {
    let two_pi: Angle = TWOPI;
    let root_two_pi: Ratio = two_pi.sqrt();
    let peak_height: PerLength = 1.0 / (sigma * root_two_pi);
    let cutoff: Length = cutoff.map_or(mm(std::f32::INFINITY), |width| width * sigma);
    move |dx: Length| -> PerLength {
        if dx.abs() < cutoff {
            let y: Ratio = dx / sigma;
            let z: Ratio = y * y;
            peak_height * (-0.5 * z).exp()
        } else {
            PerLength::ZERO
        }
    }
}

pub fn make_gauss_option(tof: Option<Tof>) -> Option<impl Fn(Length) -> PerLength> {
    tof.map(|tof| make_gauss(tof.sigma * C, Some(tof.cutoff)))
}
