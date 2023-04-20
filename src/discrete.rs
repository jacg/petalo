//! Utilities for dealing with discretization of scintillator into crystal-like
//! elements.

/// Discretization parameters of a continuous ring or tube of scintillator.
#[derive(Debug, Clone, Copy)]
pub struct Discretize {

    /// Inner radius of scintillator tube/ring
    pub r_min: Length,

    /// Thickness of scintillator tube/ring
    pub dr: Length,

    /// Axial width of scintillator elements
    pub dz: Length,

    /// Azimuthal width of scintillator elements at `r_min + dr/2`?
    pub da: Length,

    /// Place the ends of each LOR at a random position in the element
    pub smear: bool,
}

impl Discretize {

    pub fn new(r_min: Length, dr: Length, dz: Length, da: Length, smear: bool) -> Self {
        Self { r_min, dr, dz, da, smear }
    }

    /// Construct from `f32`s which are interpreted as lengths in `mm`
    pub fn from_f32s_in_mm(r_min: f32, dr: f32, dz: f32, da: f32, smear: bool) -> Self {
        Self::new(mm(r_min), mm(dr), mm(dz), mm(da), smear)
    }

    fn cell_indices(self, x: Length, y: Length, z: Length) -> (f32, f32) {
        let HelpDiscretize { d_azimuthal, .. } = self.help();
        let original_phi: Angle = y.atan2(x);
        let n_phi = ratio_(original_phi / d_azimuthal).round();
        let nz    = ratio_(z            / self.dz)    .round();
        (n_phi, nz)
    }

    /// Return a function which will adjust the position of an `xyz`-point
    /// within a scintillator element:
    /// + to the centre of the element, if `self.smear` is `false`
    /// + to a random point in the element, otherwise
    pub fn make_adjust_fn(self) -> Arc<dyn Fn(TripleLength) -> TripleLength + Send + Sync> {
        let Discretize { dr, dz, .. } = self;
        let HelpDiscretize { n_radial, d_azimuthal, .. } = self.help();

        if self.smear {
            // Move to random position in element
            let n_radial = ratio_(n_radial);
            Arc::new(move |(x, y, z)| {
                let (n_phi, n_z) = self.cell_indices(x, y, z);
                let z                   = smear(n_z)      * dz;
                let adjusted_phi: Angle = smear(n_phi)    * d_azimuthal;
                let r                   = smear(n_radial) * dr;
                let x = r * adjusted_phi.cos();
                let y = r * adjusted_phi.sin();
                (x, y, z)
            })
        } else {
            // Move to centre of element
            let r = n_radial * dr;
            Arc::new(move |(x, y, z)| {
                let (n_phi, n_z) = self.cell_indices(x, y, z);
                let z                   = n_z   * dz;
                let adjusted_phi: Angle = n_phi * d_azimuthal;
                let x = r * adjusted_phi.cos();
                let y = r * adjusted_phi.sin();
                (x, y, z)
            })
        }
    }

    // NOTE this introduces a *systematic* random bias into the sensitivity
    // image calculation. It's spread out among many elements so it might be OK,
    // but it *is* a bit fishy! Inital trials do indeed give fishy sensitivity
    // images.
    /// Generate all the points lying at the centres of the discretized elements
    /// of a detector with axial_length `scintillator_length`
    pub fn smeared_all_elements(self, scintillator_length: Length) -> impl Iterator<Item = Point> + Clone {
        let HelpDiscretize { n_radial, n_azimuthal, d_azimuthal } = self.help();
        let Discretize { dz, dr, .. } = self;

        // Ensure that there is always an element at z=0, therefore n_axial must be odd
        let half_length: Length = scintillator_length / 2.0;
        let n_half_axial: u32 = (ratio_(half_length / dz)).round() as u32;
        let n_half_axial: i32 = n_half_axial as i32;

        let     axial_indices = -n_half_axial..=n_half_axial;
        let azimuthal_indices =             0.. n_azimuthal;

        let     axial_positions = ||     axial_indices.map(move |n| n as f32 * dz);
        let azimuthal_positions = || azimuthal_indices.map(move |n| {
            let phi = d_azimuthal * n as f32;
            let r = smear(ratio_(n_radial)) * dr;
            (r * phi.cos(), r * phi.sin())
        });

        itertools::iproduct!(azimuthal_positions(), axial_positions())
            .map(|((x1,y1), z1)| Point { x: x1, y: y1, z: z1 })
    }

    /// Generate all the points lying at the centres of the discretized elements
    /// of a detector with axial_length `scintillator_length`
    pub fn centre_all_elements(self, scintillator_length: Length) -> impl Iterator<Item = Point> + Clone {
        let HelpDiscretize { n_radial, n_azimuthal, d_azimuthal } = self.help();
        let Discretize { dz, dr, .. } = self;

        let r = n_radial * dr;
        // Ensure that there is always an element at z=0, therefore n_axial must be odd
        let half_length: Length = scintillator_length / 2.0;
        let n_half_axial: u32 = (ratio_(half_length / dz)).round() as u32;
        let n_half_axial: i32 = n_half_axial as i32;

        let     axial_indices = -n_half_axial..=n_half_axial;
        let azimuthal_indices =              0..n_azimuthal;

        let     axial_positions = ||     axial_indices.map(move |n| n as f32 * dz);
        let azimuthal_positions = || azimuthal_indices.map(move |n| {
            let phi = d_azimuthal * n as f32;
            (r * phi.cos(), r * phi.sin())
        });

        itertools::iproduct!(azimuthal_positions(), axial_positions())
            .map(|((x1,y1), z1)| Point { x: x1, y: y1, z: z1 })
    }

    fn help(self) -> HelpDiscretize {
        let Discretize { r_min, dr, da, .. } = self;
        // Radius of centre of scintillator layer
        let r: Length = r_min + (dr / 2.0);
        let n_radial: Ratio = r / dr;
        let mid_circumference: Length = std::f32::consts::TAU * r;
        //let n_azimuthal: f32 = ratio_((mid_circumference / da).round::<uom::si::ratio::ratio>());
        let n_azimuthal: f32 = (ratio_(mid_circumference / da)).round();
        let d_azimuthal: Angle = TWOPI / n_azimuthal;
        HelpDiscretize {
            n_radial,
            n_azimuthal: n_azimuthal as u32,
            d_azimuthal
        }
    }
}

struct HelpDiscretize {
    /// Radial position of centre of scintillator in terms of its thickness
    n_radial: Ratio,

    /// Number of elements that fit around circumference of ring/tube
    n_azimuthal: u32,

    /// Azimuthal width of elements at `r`, after adjustment to fit
    d_azimuthal: Angle,

    // Azimithal width needs to be adjusted to fit the circumference. Axial
    // width remains the same, but we adjust the total length, if one is
    // specified. It will be specified in `all_lors_through_fov`, but not in
    // `centre_of_nearest_box_fn*`.
}

/// Randomly shift `x` to `[x - 1/2, x + 1/2)`
fn smear(x: f32) -> f32 { x - 0.5 + rand::random::<f32>() }

type TripleLength = (Length, Length, Length);
type TripleF32 = (f32, f32, f32);

/// Adapt function manipulating `uom` `Length`s into one manipulating `f32`s interpreted as `mm`
pub fn uom_mm_triplets_to_f32(f: Arc<dyn Fn(TripleLength) -> TripleLength + Sync + Send>) -> Arc<dyn Fn(TripleF32) -> TripleF32 + Sync + Send> {
    Arc::new(move |(x,y,z)| {
        let (x,y,z) = f((mm(x), mm(y), mm(z)));
        (mm_(x), mm_(y), mm_(z))
    })
}

#[cfg(test)]
mod test_discretize {
    use super::*;
    use proptest::prelude::*;
    use float_eq::assert_float_eq;
    use units::todo::Lengthf32;

    proptest! {
        #[test]
        fn generated_element_centres_invariant_when_adjusted(
            r_min in 300.0 .. (400.0 as Lengthf32),
            dr    in  10.0 .. ( 30.0 as Lengthf32),
            dz    in   1.0 .. (  5.0 as Lengthf32),
            da    in   1.0 .. (  5.0 as Lengthf32),
            l     in  30.0 .. (100.0 as Lengthf32),
        ) {
            let d = Discretize::from_f32s_in_mm(r_min, dr, dz, da, false);
            let adjust = d.make_adjust_fn();
            for Point { x, y, z } in d.centre_all_elements(mm(l)) {
                let (a,b,c) = (mm_(x), mm_(y), mm_(z));
                let (x,y,z) = adjust((x,y,z));
                let (x,y,z) = (mm_(x), mm_(y), mm_(z));
                let tol = 1e-3;
                assert_float_eq!((a,b,c), (x,y,z), abs <= (tol, tol, tol));
            }
        }
    }

    proptest!{
        #[test]
        fn smear_statistics(
            x in -30.0..30.0_f32,
        ) {
            let n = 10000;
            let sample = Vec::from_iter((0..n).map(|_| smear(x)));
            let max = *sample.iter().max_by(|a,b| a.partial_cmp(b).unwrap()).unwrap();
            let min = *sample.iter().min_by(|a,b| a.partial_cmp(b).unwrap()).unwrap();
            let sum: f32 =  sample.iter().sum();
            let spread = max - min;
            let mean = sum / n as f32;
            assert_float_eq!(spread, 1.0, abs <= 0.01);
            assert_float_eq!(mean  , x  , abs <= 0.01);
        }
    }
}

// ----- Imports ------------------------------------------------------------------------------------------
use std::sync::Arc;

use crate::Point;

use units::{
    mm, mm_, ratio_,
    Angle, Length, Ratio, TWOPI,
};
