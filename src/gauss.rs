use crate::types::{Length, TWOPI};
#[cfg(not(feature = "units"))] use crate::types::ps_to_mm;

pub use implementation::make_gauss_option;

#[cfg(not(feature = "units"))]
mod implementation {
    use super::*;

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
}

#[cfg(feature = "units")]
mod implementation {
    use super::*;

    use geometry::uom::uomcrate as guomc;
    use guomc::si::{ISQ, SI, Quantity};
    use guomc::typenum::{Z0, N1};
    use guomc::ConstZero; // num_traits::Zero;
    use geometry::uom::{Time, mm};
    use guomc::si::f32::Ratio;
    use crate::types::C;

    type PerLength = Quantity<ISQ<N1, Z0, Z0, Z0, Z0, Z0, Z0>, SI<f32>, f32>;


    // How would you make this generic over Length -> T ?
    fn make_gauss(sigma: Length, cutoff: Option<Ratio>) -> impl Fn(Length) -> PerLength {
        //let root_two_pi: Ratio = Ratio::new::<uom::si::ratio::ratio>(TWOPI.sqrt());
        let two_pi: Ratio = TWOPI;
        let root_two_pi: Ratio = two_pi.sqrt();
        let peak_height: PerLength = 1.0 / (sigma * root_two_pi);
        let cutoff: Length = cutoff.map_or(mm(std::f32::INFINITY), |width| width * sigma);
        move |dx: Length| -> PerLength {
            if dx.abs() < cutoff {
                let y: Ratio =
                    dx / sigma;
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

}
