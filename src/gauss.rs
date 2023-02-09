use units::{Angle, Length, PerLength, Ratio, TWOPI, C, mm};
use units::uom::ConstZero; // num_traits::Zero;

use crate::config::mlem::Tof;

// ---- First version: closure-based --------------------------------------------------

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

pub fn make_gauss_option_old(tof: Option<Tof>) -> Option<impl Fn(Length) -> PerLength> {
    tof.map(|tof| make_gauss(tof.sigma * C, Some(tof.cutoff)))
}

// ---- Second version: struct-based --------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct Gaussian {
    sigma: Length,
    cutoff: Length,
    peak_height: PerLength,
}

impl Gaussian {
    fn new(sigma: Length, cutoff: Option<Ratio>) -> Self {
        let two_pi: Angle = TWOPI;
        let root_two_pi: Ratio = two_pi.sqrt();
        let peak_height: PerLength = 1.0 / (sigma * root_two_pi);
        let cutoff: Length = cutoff.map_or(mm(std::f32::INFINITY), |width| width * sigma);
        Self { sigma, cutoff, peak_height }
    }

    pub fn call(&self, dx: Length) -> PerLength {
        if dx.abs() < self.cutoff {
            let y: Ratio = dx / self.sigma;
            let z: Ratio = y * y;
            self.peak_height * (-0.5 * z).exp()
        } else {
            PerLength::ZERO
        }
    }
}

pub fn make_gauss_option(tof: Option<Tof>) -> Option<Gaussian> {
    tof.map(|tof| Gaussian::new(tof.sigma * C, Some(tof.cutoff)))
}
