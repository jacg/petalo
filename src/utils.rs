use std::num::ParseFloatError;
use std::ops::{Bound, Range};

use units::{
    mm, mm_, ns, ratio, ratio_,
    Angle, Ratio, Length, TWOPI,
    todo::{Timef32, Lengthf32},
};

use crate::{Point, BoundPair, LOR};

pub fn parse_range<T: std::str::FromStr>(s: &str) -> Result<Range<T>, <T as std::str::FromStr>::Err> {
    let v = s.split("..").collect::<Vec<_>>();
    if v.len() != 2 {
        panic!("Could not find '..' when parsing range.");
    }
    let x = v[0].parse()?;
    let y = v[1].parse()?;
    Ok(x..y)
}

pub fn parse_bounds<T: std::str::FromStr>(s: &str) -> Result<BoundPair<T>, <T as std::str::FromStr>::Err> {
    let v = s.split("..").collect::<Vec<_>>();
    if v.len() != 2 {
        panic!("Could not find '..' when parsing range.");
    }

    Ok((option_to_included(parse_if_not_empty(v[0])?),
        option_to_excluded(parse_if_not_empty(v[1])?)))

}

fn option_to_included<T>(x: Option<T>) -> Bound<T> {
    if let Some(y) = x { Bound::Included(y) }
    else               { Bound::Unbounded }
}

fn option_to_excluded<T>(x: Option<T>) -> Bound<T> {
    if let Some(y) = x { Bound::Excluded(y) }
    else               { Bound::Unbounded }
}

fn parse_if_not_empty<T: std::str::FromStr>(s: &str) -> Result<Option<T>, <T as std::str::FromStr>::Err> {
    Ok(if s.is_empty() { None }
       else            { Some(s.parse()?) })
}

#[allow(clippy::many_single_char_names)]
pub fn parse_triplet<T: std::str::FromStr>(s: &str) -> Result<(T,T,T), <T as std::str::FromStr>::Err> {
    let v = s.split(',').collect::<Vec<_>>();
    assert!(v.len() == 3);
    let x = v[0].parse()?;
    let y = v[1].parse()?;
    let z = v[2].parse()?;
    Ok((x, y, z))
}

pub fn parse_lor(s: &str) -> Result<LOR, ParseFloatError> {
    let n = s.split_whitespace().collect::<Vec<_>>();
    assert!(n.len() == 8);

    let t1 = ns(n[0].parse::<Timef32>()?);
    let t2 = ns(n[1].parse::<Timef32>()?);

    let x1 = mm(n[2].parse::<Lengthf32>()?);
    let y1 = mm(n[3].parse::<Lengthf32>()?);
    let z1 = mm(n[4].parse::<Lengthf32>()?);

    let x2 = mm(n[5].parse::<Lengthf32>()?);
    let y2 = mm(n[6].parse::<Lengthf32>()?);
    let z2 = mm(n[7].parse::<Lengthf32>()?);

    let p1 = Point::new(x1, y1, z1);
    let p2 = Point::new(x2, y2, z2);
    let lor = LOR::new(t1, t2, p1, p2, ratio(1.0));
    Ok(lor)
}

// Alias to disable structopt's type magic
pub type CutoffOption<T> = Option<T>;


pub fn parse_maybe_cutoff(s: &str) -> Result<CutoffOption<Ratio>, ParseFloatError> {
    Ok(if s == "no" { None } else { Some(units::ratio(s.parse()?)) })
}

/// Group numeric digits to facilitate reading long numbers
pub fn group_digits<F: std::fmt::Display>(n: F) -> String {
    use numsep::{separate, Locale};
    separate(n, Locale::English)
}


pub mod timing {

    use super::group_digits;
    use std::time::Instant;
    use std::io::Write;

    pub struct Progress {
        previous: Instant,
    }

    impl Progress {

        #[allow(clippy::new_without_default)]
        pub fn new() -> Self { Self { previous: Instant::now() } }

        /// Print message, append ellipsis, flush stdout, stay on same line, start timer.
        pub fn start(&mut self, message: &str) {
            print!("{message} ... ");
            std::io::stdout().flush().unwrap();
            self.start_timer();
        }

        /// Print message, go to next line, start timer
        pub fn startln(&mut self, message: &str) {
            self.start(message);
            println!();
            self.start_timer();
        }

        // Print time elapsed since last start or done
        pub fn done(&mut self) {
            println!("{} ms", group_digits(self.previous.elapsed().as_millis()));
            self.start_timer();
        }

        // Print message followed by time elapsed since last start or done
        pub fn done_with_message(&mut self, message: &str) {
            println!("{message}: {} ms",
                     group_digits(self.previous.elapsed().as_millis()));
            self.start_timer();
        }

        fn start_timer(&mut self) { self.previous = Instant::now() }
    }
}


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

    /// Construct from `f32`s which are interpreted as lengths in `mm`
    pub fn from_f32s_in_mm(r_min: f32, dr: f32, dz: f32, da: f32) -> Self {
        Self { r_min: mm(r_min), dr: mm(dr), dz: mm(dz), da: mm(da) }
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

    /// Generate all the possible LORs in a detector with axial length
    /// `scintillator_length`.
    pub fn all_lors(self, scintillator_length: Length) -> impl Iterator<Item = LOR> {
        let points = self.all_element_centres(scintillator_length);
        itertools::iproduct!(points.clone(), points)
            .map(|(p1, p2)| LOR { p1, p2, dt: ns(0.0), additive_correction: ratio(1.0) })
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
        let azimuthal_indices =             0..n_azimuthal;

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
