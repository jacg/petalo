use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use petalo::weights::{Length, Time, Ratio};

// An implementation of
// https://github.com/jerenner/tofpet3d/blob/a8c3fc8293fd05547f9ca752abc259173bac57af/cc/mlem.cc#L146-L163
// in terms of petalo::weights::make_gauss
#[allow(nonstandard_style)]
fn sim_burdel(dist: Length, deltaT: Time, TOF_resolution: Time) -> Ratio {
    use petalo::weights::{make_gauss, C};
    let t = dist * 2.0 / C + deltaT ;
    let sigma = TOF_resolution;
    let cutoff = Some(3.0);
    make_gauss(sigma, cutoff)(t) * 2.0 / C
}

// Check that sim_burdel is equivalent to ale_burdel a.k.a ToFFunction
#[cfg(not(feature = "f64"))]
proptest! {
    #[test]
    fn gauss_equivalence(
        dist       in -100.0 .. (100.0 as Length),
        dt         in -100.0 .. (100.0 as Time),
        resolution in   20.0 .. (200.0 as Time),
    ) {
        use cmlem::ale_burdel;
        let sane = sim_burdel(dist, dt, resolution);
        let mad  = ale_burdel(dist, dt, resolution);
        assert_approx_eq!(sane, mad);
    }
}
