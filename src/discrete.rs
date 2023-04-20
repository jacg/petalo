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
}

impl Discretize {

    pub fn new(r_min: Length, dr: Length, dz: Length, da: Length) -> Self {
        Self { r_min, dr, dz, da }
    }

    /// Construct from `f32`s which are interpreted as lengths in `mm`
    pub fn from_f32s_in_mm(r_min: f32, dr: f32, dz: f32, da: f32) -> Self {
        Self::new(mm(r_min), mm(dr), mm(dz), mm(da))
    }

    /// Return a function which will move an `xyz`-point to the centre of the
    /// nearest discrete element. Inputs and outputs are `uom` `Length`s. See
    /// also `centre_of_nearest_box_fn_f32`.
    pub fn centre_of_nearest_box_fn_uom(self) -> impl Fn(Length, Length, Length) -> (Length, Length, Length) + Copy {
        let HelpDiscretize { r, d_azimuthal, .. } = self.help();
        move |x:Length, y:Length, z:Length| {
            let original_phi: Angle = y.atan2(x);
            let n = ratio_(original_phi / d_azimuthal).round();
            let adjusted_phi: Angle = n * d_azimuthal;
            let x = r * adjusted_phi.cos();
            let y = r * adjusted_phi.sin();
            let z = (ratio_(z / self.dz)).round() * self.dz;
            (x, y, z)
        }
    }

    /// Return a function which will move an `xyz`-point to the centre of the
    /// nearest discrete element. Inputs and outputs are `f32`s interpreted as
    /// lengths in `mm`. See also `centre_of_nearest_box_fn_uom`.
    pub fn centre_of_nearest_box_fn_f32(self) -> impl Fn(f32, f32, f32) -> (f32, f32, f32) + Copy {
        let adjust = self.centre_of_nearest_box_fn_uom();
        move |x,y,z| {
            let (x,y,z) = adjust(mm(x), mm(y), mm(z));
            (mm_(x), mm_(y), mm_(z))
        }
    }

    /// Generate all the possible pairs of `Point`s at the centres of the
    /// scintillator elements in a detector with axial length `scintillator_length`.
    /// NB: the number of LORs quickly becomes intractable with increasing detector size.
    pub fn all_lors(self, scintillator_length: Length) -> impl Iterator<Item = (Point, Point)>
    where
    {
        let points = self.all_element_centres(scintillator_length);
        itertools::iproduct!(points.clone(), points)
            .filter(|(p1,p2)| p1 != p2)
    }

    /// Generate all the points lying at the centres of the discretized elements
    /// of a detector with axial_length `scintillator_length`
    pub fn all_element_centres(self, scintillator_length: Length) -> impl Iterator<Item = Point> + Clone {
        let HelpDiscretize { r, n_azimuthal, d_azimuthal } = self.help();
        let Discretize { dz, .. } = self;

        // Ensure that there is always an element at z=0, therefore n_axial must be odd
        let half_length: Length = scintillator_length / 2.0;
        let n_half_axial: u32 = (ratio_(half_length / dz)).round() as u32;
        let n_half_axial: i32 = n_half_axial as i32;

        let     axial_indices = -n_half_axial..=n_half_axial;
        let azimuthal_indices =             0.. n_azimuthal;

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
        let mid_circumference: Length = std::f32::consts::TAU * r;
        //let n_azimuthal: f32 = ratio_((mid_circumference / da).round::<uom::si::ratio::ratio>());
        let n_azimuthal: f32 = (ratio_(mid_circumference / da)).round();
        let d_azimuthal: Angle = TWOPI / n_azimuthal;
        HelpDiscretize {
            r,
            n_azimuthal: n_azimuthal as u32,
            d_azimuthal
        }
    }
}

struct HelpDiscretize {
    /// Radial position of centre of scintillator
    r: Length,

    /// Number of elements that fit around circumference of ring/tube
    n_azimuthal: u32,

    /// Azimuthal width of elements at `r`, after adjustment to fit
    d_azimuthal: Angle,

    // Azimithal width needs to be adjusted to fit the circumference. Axial
    // width remains the same, but we adjust the total length, if one is
    // specified. It will be specified in `all_lors_through_fov`, but not in
    // `centre_of_nearest_box_fn*`.
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
            let d = Discretize::from_f32s_in_mm(r_min, dr, dz, da);
            let adjust = d.centre_of_nearest_box_fn_uom();
            for Point { x, y, z } in d.all_element_centres(mm(l)) {
                let (a,b,c) = (mm_(x), mm_(y), mm_(z));
                let (x,y,z) = adjust(x,y,z);
                let (x,y,z) = (mm_(x), mm_(y), mm_(z));
                let tol = 1e-3;
                assert_float_eq!((a,b,c), (x,y,z), abs <= (tol, tol, tol));
            }
        }
    }
}

// ----- Imports ------------------------------------------------------------------------------------------
use crate::Point;

use units::{
    mm, mm_, ratio_,
    Angle, Length, TWOPI,
};
