// NOTE only using first vertex, for now
#[allow(nonstandard_style)]
pub (crate) fn lor_from_discretized_vertices(d: &Reco) -> impl Fn(&[Vertex]) -> Option<Hdf5Lor> + Send + Sync {
    let &Reco::Discrete { r_min, dr, dz, da, smear } = d else {
        panic!("lor_from_discretized_vertices called with variant other than Reco::Discrete")
    };
    use petalo::discrete::{Discretize, uom_mm_triplets_to_f32};
    let adjust = uom_mm_triplets_to_f32(Discretize { r_min, dr, dz, da, smear }.make_adjust_fn());
    move |vertices| {
        let mut in_scint = vertices_in_scintillator(vertices);

        // Optimistic version: grabs all remaining KE in first vertex
        let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: E2, ..} = in_scint.find(|v| v.track_id == 2)?;
        let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: E1, ..} = in_scint.find(|v| v.track_id == 1)?;

        // // Pessimistic version: picks vertex with highest dKE
        // use ordered_float::NotNan;
        // let in_scint = &mut vertices_in_scintillator(vertices); // Working around filter consuming the iterator
        // let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: b2, post_KE: a2, ..} = in_scint.filter(|v| v.track_id == 2).max_by_key(|v| NotNan::new(v.pre_KE - v.post_KE).unwrap())?;
        // let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: b1, post_KE: a1, ..} = in_scint.filter(|v| v.track_id == 1).max_by_key(|v| NotNan::new(v.pre_KE - v.post_KE).unwrap())?;
        // let (E1, E2) = (b1 - a1, b2 - a2);

        // let (ox1, oy1, oz1, ox2, oy2, oz2) = (x1, y1, z1, x2, y2, z2);
        let (x1, y1, z1) = adjust((x1, y1, z1));
        let (x2, y2, z2) = adjust((x2, y2, z2));

        // println!();
        // let dp1 = dist((ox1, oy1, oz1), (x1, y1, z1));
        // let dp2 = dist((ox2, oy2, oz2), (x2, y2, z2));
        // let mean_z = (z2+z1) / 2.0;
        // if E1.min(E2) < 500.0 { return None }
        // println!("{ox1:6.1} {oy1:6.1} {oz1:6.1}   {ox2:6.1} {oy2:6.1} {oz2:6.1}     {dp1:4.1}  {dp2:4.1}   {mean_z:6.1}");
        // println!("{x1:6.1} {y1:6.1} {z1:6.1}   {x2:6.1} {y2:6.1} {z2:6.1}    {E1:5.1} {E2:5.1}");

        Some(Hdf5Lor { dt: t2 - t1, x1, y1, z1, x2, y2, z2, q1: f32::NAN, q2: f32::NAN, E1, E2 })
    }
}

// ----- Imports -----------------------------------------------------------------------------------------
use crate::{Reco, vertices_in_scintillator};
use petalo::io::hdf5::{Hdf5Lor, mc::Vertex};

// ----- TESTS ------------------------------------------------------------------------------------------
#[cfg(test)]
mod test_discretize {
    use rstest::rstest;
    use float_eq::assert_float_eq;

    use std::f32::consts::SQRT_2 as ROOT2;
    use std::f32::consts::TAU;

    #[rstest]
    //              rmin  dr     dz    da       x-in    y-in   z-in     x-out   y-out  z-out
    // On x/y-axis
    #[case::x_pos(( 90.0, 20.0,  3.0,  TAU), ( 109.9,    0.2,  28.6), ( 100.0 ,   0.0, 30.0))]
    #[case::x_neg(( 90.0, 20.0,  4.0,  TAU), (- 91.3,    0.2,  33.9), (-100.0 ,   0.0, 32.0))]
    #[case::y_pos(( 90.0, 20.0,  5.0,  TAU), (   0.2,   95.3,  50.1), (   0.0 , 100.0, 50.0))]
    #[case::y_neg(( 90.0, 20.0,  5.5,  TAU), (   0.2, -109.9,  52.3), (   0.0 ,-100.0, 55.0))]
    // Check that quadrants are preserved
    #[case::quad1((  1.9,  0.2,  1.0, 0.01), ( ROOT2,  ROOT2, -23.4), ( ROOT2, ROOT2, -23.0))]
    #[case::quad2((  1.9,  0.2,  1.0, 0.01), (-ROOT2,  ROOT2, - 9.9), (-ROOT2, ROOT2, -10.0))]
    #[case::quad3((  1.9,  0.2,  1.0, 0.01), (-ROOT2, -ROOT2, - 9.9), (-ROOT2,-ROOT2, -10.0))]
    #[case::quad4((  1.9,  0.2,  1.0, 0.01), ( ROOT2, -ROOT2,  23.4), ( ROOT2,-ROOT2,  23.0))]
    // Near x-axis
    #[case::xish1((280.0, 30.0,  4.0, 2.0 ), ( 309.9,    0.0, -31.6), ( 295.0,   0.0, -32.0))]
    #[case::xish2((280.0, 30.0,  5.0, 2.0 ), ( 280.1,    2.8,  31.6), ( 295.0,   2.0,  30.0))]
    #[case::xish3((280.0, 30.0,  6.0, 1.5 ), ( 294.6,   -0.8, -31.6), ( 295.0,  -1.5, -30.0))]
    // Other cases are very fiddly to verify by hand, but we should add some, in principle
    fn test_nearest_centre_of_box(
        #[case] detector: (f32, f32, f32, f32),
        #[case] vertex  : (f32, f32, f32),
        #[case] expected: (f32, f32, f32),
    ) {
        use petalo::discrete::{Discretize, uom_mm_triplets_to_f32};
        let (rmin, dr, dz, da) = detector;
        let move_to_centre_of_element = uom_mm_triplets_to_f32(Discretize::from_f32s_in_mm(rmin, dr, dz, da, false).make_adjust_fn());
        let (x, y, z) = move_to_centre_of_element(vertex);
        assert_float_eq!((x,y,z), expected, abs <= (0.01, 0.01, 0.01));
    }
}
