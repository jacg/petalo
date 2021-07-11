use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use petalo::types::{Length, Time, Ratio};

// An implementation of
// https://github.com/jerenner/tofpet3d/blob/a8c3fc8293fd05547f9ca752abc259173bac57af/cc/mlem.cc#L146-L163
// in terms of petalo::weights::make_gauss
#[allow(nonstandard_style)]
fn sim_burdel(dist: Length, dt: Time, sigma: Time) -> Ratio {
    use petalo::types::C;
    use petalo::gauss::make_gauss_option;
    let x = dist + dt * C / 2.0;
    let cutoff = Some(3.0);
    let gauss = make_gauss_option(Some(sigma / 2.0), cutoff).unwrap();
    gauss(x)
}

// Check that sim_burdel is equivalent to ale_burdel a.k.a ToFFunction
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
        println!("ratio: {}", sane / mad);
        assert_approx_eq!(sane, mad);
    }
}
