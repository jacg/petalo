#[allow(nonstandard_style)]
pub (crate) fn lor_from_discretized_vertices(d: &Reco) -> impl Fn(&[Vertex]) -> Option<Hdf5Lor> + Send + Sync {
    let &Reco::Discrete { r_min, dr, dz, da, adjust } = d else {
        panic!("lor_from_discretized_vertices called with variant other than Reco::Discrete")
    };
    //let adjust = uom_mm_triplets_to_f32(Discretize { r_min, dr, dz, da, smear }.make_adjust_fn());
    let discretize = Discretize { r_min, dr, dz, da, adjust };
    move |vertices| {
        let in_scint = vertices_in_scintillator(vertices);

        // // Optimistic version: grabs all remaining KE in first vertex
        // let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: E2, ..} = in_scint.find(|v| v.track_id == 2)?;
        // let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: E1, ..} = in_scint.find(|v| v.track_id == 1)?;

        // // Pessimistic version: picks vertex with highest dKE
        // use ordered_float::NotNan;
        // let in_scint = &mut vertices_in_scintillator(vertices); // Working around filter consuming the iterator
        // let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: b2, post_KE: a2, ..} = in_scint.filter(|v| v.track_id == 2).max_by_key(|v| NotNan::new(v.pre_KE - v.post_KE).unwrap())?;
        // let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: b1, post_KE: a1, ..} = in_scint.filter(|v| v.track_id == 1).max_by_key(|v| NotNan::new(v.pre_KE - v.post_KE).unwrap())?;
        // let (E1, E2) = (b1 - a1, b2 - a2);

        // // Smear positions
        // let (x1, y1, z1) = adjust((x1, y1, z1));
        // let (x2, y2, z2) = adjust((x2, y2, z2));

        // Some(Hdf5Lor { dt: t2 - t1, x1, y1, z1, x2, y2, z2, q1: f32::NAN, q2: f32::NAN, E1, E2 })

        let (vs1, vs2): (Vec<_>, _) = in_scint
            .filter   (|v| v.track_id <  3)
            .partition(|v| v.track_id == 1);

        let (E1, (x1, y1, z1)) = box_with_higest_total_energy_centre_doi(&vs1, discretize)?;
        let (E2, (x2, y2, z2)) = box_with_higest_total_energy_centre_doi(&vs2, discretize)?;

        Some(Hdf5Lor { dt: 0.0, x1, y1, z1, x2, y2, z2, q1: f32::NAN, q2: f32::NAN, E1, E2 })

    }
}

/// Find the box with the highest deposited energy.
///
/// Return:
/// + total energy deposited in that box
/// + `(x, y, z)` coordinate of
///   - centre of the box in `z` and `phi`
///   - barycentre in `r` of vertices in box, weighted by deposited energy
fn box_with_higest_total_energy_centre_doi(vertices: &[Vertex], discretize: Discretize) -> Option<(f32, (f32, f32, f32))> {
    use itertools::Itertools;
    let zero_energy = keV(0.0);
    let zero_energy_length = mm(0.0) * zero_energy;
    vertices.iter()
        .cloned()
        .map(|Vertex { x, y, z, pre_KE, post_KE, .. }| ((mm(x),mm(y),mm(z)), keV(pre_KE - post_KE)))
        .into_grouping_map_by(|&((x,y,z), _energy)| discretize.cell_indices(x, y, z))
        .fold(
            (zero_energy, zero_energy_length),   // Initialize accumulator
            |(e_acc, re_acc), _i, ((x,y,_), e)|
            ( e_acc + e             ,  // accumulate energy
             re_acc + e * x.hypot(y))  // accumulate energy-weighted r
        )
        .into_iter()
        .max_by(|(_,(e1,_)), (_,(e2,_))| e1.partial_cmp(e2).unwrap())
        .map(|(i, (e, re))| {
            let (x,y,z) = discretize.indices_to_centre(i);
            let r_centre: Length = x.hypot(y);
            let r_bary: Length = re / e;
            let scale = r_bary / r_centre;
            let scale = units::ratio(1.0);
            (keV_(keV(511.0)), (mm_(x*scale),
                       mm_(y*scale),
                       mm_(z      )))
        })
}


// ----- Imports -----------------------------------------------------------------------------------------
use crate::{Reco, vertices_in_scintillator};
use petalo::{
    discrete::Discretize,
    io::hdf5::{Hdf5Lor, mc::Vertex},
};
use units::{mm, mm_, keV, keV_, Length};

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
        let move_to_centre_of_element = uom_mm_triplets_to_f32(
            Discretize::from_f32s_in_mm(rmin, dr, dz, da, petalo::discrete::Adjust::Centre)
                .make_adjust_fn());
        let (x, y, z) = move_to_centre_of_element(vertex);
        assert_float_eq!((x,y,z), expected, abs <= (0.01, 0.01, 0.01));
    }
}
