//! Configuration file parser for MLEM

use std::{fs, fmt::{Display, Debug}};
use std::str::FromStr;
use std::path::PathBuf;

use serde::{Deserialize, Deserializer, de};


use crate::{Length, Ratio, Time};

#[cfg(test)]
fn deserialize_uom_opt<'d, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'d>,
    T: FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    Option::<&str>::deserialize(deserializer)?
        .map(str::parse::<T>)
        .transpose()
        .map_err(de::Error::custom)
}

fn deserialize_uom<'d, D, T>(deserializer: D) -> Result<T, D::Error>
where
    D: Deserializer<'d>,
    T: FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    Option::<&str>::deserialize(deserializer)?
        .map(str::parse::<T>)
        .unwrap()
        .map_err(de::Error::custom)
}

fn _deserialize_uom_3d_opt<'d, D, T>(deserializer: D) -> Result<Option<(T, T, T)>, D::Error>
where
    D: Deserializer<'d>,
    T: FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    Option::<(&str, &str, &str)>::deserialize(deserializer)?
        .map(|(x,y,z)| tr_tup_res((x.parse(), y.parse(), z.parse())))
        .transpose()
        .map_err(de::Error::custom)
}

fn deserialize_uom_3d<'d, D, T>(deserializer: D) -> Result<(T, T, T), D::Error>
where
    D: Deserializer<'d>,
    T: FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    Option::<(&str, &str, &str)>::deserialize(deserializer)?
        .map(|(x,y,z)| tr_tup_res((x.parse(), y.parse(), z.parse())))
        .unwrap()
        .map_err(de::Error::custom)
}


/// Transpose 3-tuple of `Result`
///
/// `Ok` if all elements `Ok`; if any element is an `Err` return the first one.
///
/// # Examples
/// `(Ok(a),  Ok(b),  Ok(c)) -> Ok((a, b, c))`
/// `(Ok(a), Err(b),  Ok(c)) -> Err(b)`
/// `(Ok(a), Err(b), Err(c)) -> Err(b)`
fn tr_tup_res<O, E>((x,y,z): (Result<O, E>, Result<O, E>, Result<O, E>)) -> Result<(O, O, O), E> {
    Ok((x?, y?, z?))
}

#[derive(Deserialize, Debug, Default)]
#[serde(deny_unknown_fields)]
pub struct Config {

    /// Specification of data to be used to reconstruct image
    #[serde(default = "mandatory")]
    pub input: Input,

    /// Number of MLEM/OSEM iterations/subsets to use in reconstruction
    #[serde(default = "mandatory")]
    pub iterations: Iterations,

    /// Field of view (FOV) size: number of voxels and physical dimensions
    #[serde(default = "mandatory")]
    pub fov: Fov,

    /// Time of Flight (TOF) parameters
    pub tof: Option<Tof>,

    /// Sensitivity image to use in scatter correction
    pub attenuation_correction: Option<AttenuationCorrection>,

    pub scatter_correction: Option<Scatter>,
}

#[derive(Deserialize, Debug, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct Input {

    /// HDF5 file containing reconstructed LORs from which to reconstruct image
    #[serde(default = "mandatory")]
    pub file: PathBuf,

    /// The dataset location inside the input file
    #[serde(default = "mandatory")]
    pub dataset: String,

    #[serde(default)]
    pub energy: Bounds<f32>,

    #[serde(default)]
    pub charge: Bounds<f32>,

    #[serde(default)]
    pub events: Bounds<usize>,

}

#[derive(Deserialize, Debug, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct Bounds<T> {

    pub min: Option<T>, // TODO use uom units
    pub max: Option<T>,

}

impl<T> Bounds<T> {
    pub fn new(min: Option<T>, max: Option<T>) -> Self { Self { min, max } }
    pub fn none() -> Self { Self::new(None, None) }
}

impl<T: PartialOrd + Copy> Bounds<T> {
    pub fn contains(&self, e: T) -> bool {
        (self.min.is_none() || self.min.unwrap() <= e) &&
        (self.max.is_none() || self.max.unwrap() >  e)
    }
}

#[derive(Deserialize, Debug, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct Iterations {

    /// Number of MLEM or OSEM iterations to perform
    #[serde(default = "mandatory")]
    pub number: usize,

    /// Number of OSEM subsets to use. For MLEM, `subsets` = 1
    #[serde(default = "mandatory")]
    pub subsets: usize,

}

#[derive(Deserialize, Debug, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct Fov {

    #[serde(default = "mandatory")]
    pub nvoxels: (usize, usize, usize),

    #[serde(default = "mandatory")]
    #[serde(deserialize_with = "deserialize_uom_3d")]
    pub size: (Length, Length, Length),

}

#[derive(Deserialize, Debug, Clone, Copy, Default)]
#[serde(deny_unknown_fields)]
pub struct Tof {

    #[serde(default)]
    #[serde(deserialize_with = "deserialize_uom")]
    pub sigma: Time,

    #[serde(default = "three")]
    pub cutoff: Ratio,

}

fn three() -> Ratio { geometry::units::ratio(3.0) }

#[derive(Deserialize, Debug, Clone, Default)]
#[serde(deny_unknown_fields)]
pub struct AttenuationCorrection {

    #[serde(default)]
    pub sensitivity_image: PathBuf,

}

#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Scatter {
    pub phi: Option<Bins>,
    pub   r: Option<BinsMax<Length>>,
    pub  dz: Option<BinsMax<Length>>,
    pub  dt: Option<BinsMax<Time>>,
    pub   z: Option<BinsLength>,
}

#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Bins {
    pub bins: usize,
}

#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct BinsMax<T>
where
    T : FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    pub bins: usize,
    #[serde(deserialize_with = "deserialize_uom")]
    pub max: T,
}

#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct BinsLength {
    pub bins: usize,
    #[serde(deserialize_with = "deserialize_uom")]
    pub length: Length,
}

pub fn read_config_file(path: PathBuf) -> Config {
    let config: String = fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("Couldn't read config file `{:?}`", path));
    toml::from_str(&config).unwrap()
}



// Hack to allow mandatory fields to be missing during testing.
#[cfg(not(test))]
fn mandatory<T>() -> T { panic!("MISSING MANDATORY FIELD. TODO: which one?") }
#[cfg(test)]
fn mandatory<T: Default>() -> T { T::default() }



#[cfg(test)]
mod tests {
    use super::*;

    use geometry::units::{cm, mm, ps, ratio};

    // ----- Test an example on-disk config file -----------------------------------------
    #[test]
    fn test_config_file() {
        let config = read_config_file("test-data/mlem-config.toml".into());

        assert_eq!(config.input.file   , PathBuf::from_str("data/some-lors.h5").unwrap());
        assert_eq!(config.input.dataset, String ::from    ("reco_info/lors"));

        assert_eq!(config.iterations.number ,  4);
        assert_eq!(config.iterations.subsets, 20);

        let tof = config.tof.unwrap();
        assert_eq!(tof.sigma, ps(200.0));
        assert_eq!(tof.cutoff, ratio(3.0));

        let scatter = config.scatter_correction.unwrap();
        let phi = scatter.phi.unwrap();
        let r   = scatter.r  .unwrap();
        let dz  = scatter.dz .unwrap();
        let dt  = scatter.dt .unwrap();
        let  z  = scatter. z .unwrap();

        assert_eq!(phi.bins,    30   );

        assert_eq!(  r.bins,    31   );
        assert_eq!(  r.max , mm(32.0));

        assert_eq!( dz.bins,    33   );
        assert_eq!( dz.max , mm(34.0));

        assert_eq!( dt.bins,    35   );
        assert_eq!( dt.max , ps(36.0));

        assert_eq!(  z.bins,      37   );
        assert_eq!(  z.length, cm(38.0));
    }

    // ----- Some helpers to make the tests more concise ---------------------------------
    //  ---  Parse string as TOML  -------------------------
    fn parse<'d, D: Deserialize<'d>>(input: &'d str) -> D {
        toml::from_str(input).unwrap()
    }
    // //  ---  Parse string as TOML, with explicit error reporting -------------------------
    fn parse_carefully<'d, D: Deserialize<'d>>(input: &'d str) -> Result<D, toml::de::Error> {
        toml::from_str(input)
    }
    fn parse_config<'d>(input: &'d str) -> Result<Config, toml::de::Error> {
        parse_carefully(input)
    }
    //  ---  Macro for concise assertions about vlues of parsed fields -------------------
    macro_rules! check {
        ($type:ident($text:expr).$field:ident = $expected:expr) => {
            let config: $type = parse::<$type>($text);
            println!("DESERIALIZED: {config:?}");
            assert_eq!(config.$field, $expected);
        };
        ($type:ident($text:expr) fields: $($field:ident = $expected:expr);+$(;)?) => {
            let config: $type = parse::<$type>($text);
            println!("DESERIALIZED: {config:?}");
            $(assert_eq!(config.$field, $expected);)*
        }
    }
    // ----- Test deserializing of individual aspects of the Config type ----------------
    // ----- Make sure that unknown fields are not accepted -----------------------------
    #[test]
    #[should_panic]
    fn config_reject_unknown_field() {
        parse("unknown_field = 666")
    }
    // ----- LORs -----------------------------------------------------------------------
    #[test]
    fn config_input_file() {
        let input = parse::<Config>(r#"
            [input]
            file    = "some/file.h5"
            dataset = "some/dataset"
        "#).input;
        assert_eq!(input.file   , PathBuf::from_str("some/file.h5").unwrap());
        assert_eq!(input.dataset, String ::from    ("some/dataset"));
        assert_eq!(input.energy.min, None);
        assert_eq!(input.energy.max, None);
        assert_eq!(input.charge.min, None);
        assert_eq!(input.charge.max, None);
        assert_eq!(input.events.min, None);
        assert_eq!(input.events.max, None);
    }

    #[test]
    fn config_input_file_energy_cut() {
        let input = parse::<Config>(r#"
            [input]
            energy.min = 123
            charge.min = 444
            events.min = 100_000
        "#).input;
        assert_eq!(input.energy.min, Some(123.0));
        assert_eq!(input.energy.max, None);
        assert_eq!(input.charge.min, Some(444.0));
        assert_eq!(input.charge.max, None);
        assert_eq!(input.events.min, Some(100_000));
        assert_eq!(input.events.max, None);

        let input = parse::<Config>(r#"
            [input]
            energy.max = 456
            charge.max = 999
            events.max = 123
        "#).input;
        assert_eq!(input.energy.min, None);
        assert_eq!(input.energy.max, Some(456.0));
        assert_eq!(input.charge.min, None);
        assert_eq!(input.charge.max, Some(999.0));
        assert_eq!(input.events.min, None);
        assert_eq!(input.events.max, Some(123));

        let input = parse::<Config>(r#"
            [input]
            energy.min = 434
            energy.max = 600
            charge.min = 111
            charge.max = 222
            events.min = 100_000_000
            events.max = 200_000_000
        "#).input;
        assert_eq!(input.energy.min, Some(434.0));
        assert_eq!(input.energy.max, Some(600.0));
        assert_eq!(input.charge.min, Some(111.0));
        assert_eq!(input.charge.max, Some(222.0));
        assert_eq!(input.events.min, Some(100_000_000));
        assert_eq!(input.events.max, Some(200_000_000));
    }
    // ----- iterations -----------------------------------------------------------------
    #[test]
    fn config_iterations() {
        let iterations = parse::<Config>(r#"
             [iterations]
             number = 50
             subsets = 1
        "#).iterations;

        assert_eq!(iterations.number, 50);
        assert_eq!(iterations.subsets, 1);

        let iterations = parse::<Config>(r#"iterations = { number = 4, subsets = 20 }"#).iterations;

        assert_eq!(iterations.number,   4);
        assert_eq!(iterations.subsets, 20);
    }
    // ----- Test FOV parameters ---------------------------------------------------------
    #[test]
    fn config_fov() {
        let fov = parse::<Config>(r#"
                     [fov]
                     nvoxels = [  10    ,   20    ,  30    ]
                     size    = ["123 mm", "456 mm", "78 cm"]
               "#).fov;
        assert_eq!(fov.nvoxels, (    10   ,    20    ,    30   ));
        assert_eq!(fov.size   , (mm(123.0), mm(456.0), cm(78.0)));

    }
    // ----- Test TOF parameters ---------------------------------------------------------
    #[test]
    fn config_tof() {
        let tof = parse::<Config>(r#"
                     [tof]
                     sigma = "200 ps"
                     cutoff = 3
               "#).tof.unwrap();
        assert_eq!(tof.sigma , ps(200.0));
        assert_eq!(tof.cutoff, ratio(3.0));
    }

    #[test]
    fn config_tof_missing() {
        let tof = parse::<Config>("").tof;
        assert!(tof.is_none());
    }

    #[test]
    fn config_tof_cutoff_default() {
        let tof = parse::<Config>(r#"
                     [tof]
                     sigma = "123.4 ps"
               "#).tof.unwrap();
        assert_eq!(tof.sigma , ps(123.4));
        assert_eq!(tof.cutoff, ratio(3.0));
    }
    // ----- Test Scattergram parameters ------------------------------------------------
    #[test]
    fn config_attenuation_correction_yes() {
        let corr = parse::<Config>(r#"
                      [attenuation_correction]
                      sensitivity_image = "some/sensitivity_image.raw"
               "#).attenuation_correction.unwrap();
        assert_eq!(corr.sensitivity_image, PathBuf::from_str("some/sensitivity_image.raw").unwrap());
    }

    #[test]
    fn config_attenuation_correction_no() {
        let corr = parse::<Config>("").attenuation_correction;
        assert!(corr.is_none());
    }
    // ----- Test attenuation correction parameters -------------------------------------
    #[test]
    fn config_scattergram() {
        let c: Config = parse(r#"
                 [scatter_correction]
                 phi.bins = 12
                 r  .bins = 34
                 dz .bins = 99
                 dt .bins = 98
                 z  .bins = 97

                 r  .max  = "56 mm"
                 dz .max  = "78 mm"
                 dt .max  = "90 ps"
                 z.length = "38 cm"
              "#);
        let scatter = c.scatter_correction.unwrap();

        let phi = scatter.phi.unwrap();
        let r   = scatter.r  .unwrap();
        let dz  = scatter.dz .unwrap();
        let dt  = scatter.dt .unwrap();
        let  z  = scatter. z .unwrap();

        assert_eq!(phi.bins  ,    12   );

        assert_eq!(  r.bins  ,    34   );
        assert_eq!(  r.max   , mm(56.0));

        assert_eq!( dz.bins  ,    99   );
        assert_eq!( dz.max   , mm(78.0));

        assert_eq!( dt.bins  ,    98   );
        assert_eq!( dt.max   , ps(90.0));

        assert_eq!(  z.bins  ,    97   );
        assert_eq!(  z.length, cm(38.0));
    }
    // -----------------------------------------------------------------------------------
    // The tests that follow should be read in order: they tell the story of why
    // and how we need to jump through a number of hoops in order to parse uom
    // values.

    // uom types are, by default, deserialized from TOML numbers, with the
    // Quantity's base unit being inferred!
    #[test]
    fn toml_uom_without_units() {
        #[derive(Deserialize, Debug)]
        struct X { a: Option<Time> }
        check!(X("a = 2").a = Some(ps(2.))); // NOTE: no `ps` in input
        check!(X("     ").a = None        );
    }
    // This defeats the point of uom: making units explicit.
    // So we need to work a bit harder to have the units be parsed ...

    // By default, serde treats Option fields as ones that may be missing
    #[test]
    fn toml_u32() {
        #[derive(Deserialize, Debug)]
        struct X { a: Option<u32> }
        check!(X("a = 10").a = Some(10));
        check!(X("      ").a = None    );
    }

    // But, if using `deserialize_with` the rules change, and missing fields
    // generate an error.
    #[test]
    #[should_panic] // TODO: replace unwrap with match
    fn toml_u32_deserialize_with_missing_field() {
        #[derive(Deserialize, Debug)]
        struct X {
            #[serde(deserialize_with = "bbb")]
            a: Option<u32>
        }
        check!(X("      ").a = None);

        fn bbb<'d, D: Deserializer<'d>>(deserializer: D) -> Result<Option<u32>, D::Error> {
            Ok(Option::<u32>::deserialize(deserializer)?)
        }
    }
    // See
    //
    // https://stackoverflow.com/questions/44301748/how-can-i-deserialize-an-optional-field-with-custom-functions-using-serde
    //
    // Consequetly, me must add `[serde(default)]` for the field to be truly optional
    #[test]
    fn toml_u32_deserialize_with() {
        #[derive(Deserialize, Debug)]
        struct X {
            #[serde(default)]
            #[serde(deserialize_with = "bbb")]
            a: Option<u32>
        }
        check!(X("a = 10").a = Some(10));
        check!(X("      ").a = None    );

        fn bbb<'d, D: Deserializer<'d>>(deserializer: D) -> Result<Option<u32>, D::Error> {
            Ok(Option::<u32>::deserialize(deserializer)?)
        }
    }

    // The other problem with parsing uom-quantities-with-units from TOML, is
    // related to TOML understanding very few types: numbers, strings, maps,
    // arrays, dates and not much more.
    //
    // Something like `2 ps` cannot be parsed as a number, so we must stick
    // quotes around it in the TOML source, and first parse it as a string, then
    // parse that string into the relevant `uom` type
    #[test]
    fn toml_time_with_units_via_str() {
        #[derive(Deserialize, Debug)]
        struct X {
            #[serde(default)]
            #[serde(deserialize_with = "parse_uom_time")]
            a: Option<Time>
        }
        check!(X(r#"a = "2 ps""#).a = Some(ps(2.   )));
        check!(X(r#"a = "2 ns""#).a = Some(ps(2000.)));
        check!(X(r#"          "#).a = None           );

        fn parse_uom_time<'d, D: Deserializer<'d>>(deserializer: D) -> Result<Option<Time>, D::Error> {
            Option::<&str>::deserialize(deserializer)?
                .map(str::parse::<Time>)
                .transpose()
                .map_err(de::Error::custom)
        }
    }

    // ... which relies on the parsers provided by `uom`
    #[test]
    fn uom_parse() -> Result<(), Box<dyn std::error::Error>> {
        let t: Time = "2 ps".parse()?;
        assert_eq!(t, ps(2.));

        let t: Time = "2 ns".parse()?;
        assert_eq!(t, ps(2000.));
        Ok(())
    }

    // `parse_uom_time` used above, works only for `Time`. We need a generic
    // version, which is implemented outside of the test module and tested here
    #[test]
    fn toml_with_units_generic_deserialize() {
        #[derive(Deserialize, Debug)]
        struct X {
            #[serde(default)]
            #[serde(deserialize_with = "deserialize_uom_opt")]
            t: Option<Time>,
            #[serde(default)]
            #[serde(deserialize_with = "deserialize_uom_opt")]
            l: Option<Length>,
        }
        check!(X(r#"t = "2 ps""#).t = Some(ps(2.)));
        check!(X(r#"l = "3 mm""#).l = Some(mm(3.)));
        check!(X(r#"          "#).t = None        );
    }

}


// TODO: consider making the representation readable by the config file parser
impl Display for Config {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("\n[input]\n{}"                  , self.input))?;
        f.write_fmt(format_args!("\n[iterations]\n{}"             , self.iterations))?;
        f.write_fmt(format_args!("\n[fov]\n{}"                    , self.fov))?;

        f.write_str("\n\n[tof]\n")?;
        if let Some(Tof { sigma, cutoff }) = &self.tof {
            f.write_fmt(format_args!("sigma = {sigma:?}\n"))?;
            f.write_fmt(format_args!("cutoff = {cutoff:?}"))?;
        } else {
            f.write_str("OFF")?;
        }

        f.write_str("\n\n[attenuation_correction]\n")?;
        if let Some(AttenuationCorrection { sensitivity_image }) = &self.attenuation_correction {
            f.write_fmt(format_args!("{}" , sensitivity_image.display()))?;
        } else {
            f.write_str("OFF")?;
        }

        f.write_str("\n\n[scatter_correction]\n")?;
        if let Some(scatter) = &self.scatter_correction {
            f.write_fmt(format_args!("{}" , scatter))?;
        } else {
            f.write_str("OFF")?;
        }
        f.write_str("\n")
    }
}

impl Display for Input {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("file    = {}\n", self.file.display()))?;
        f.write_fmt(format_args!("dataset = {}\n", self.dataset))?;
        f.write_fmt(format_args!("energy  : {}\n", self.energy))?;
        f.write_fmt(format_args!("charge  : {}\n", self.charge))?;
        f.write_fmt(format_args!("events  : {}\n", self.events))?;
        Ok(())
    }
}

#[allow(clippy::useless_format)]
impl<T: Display + Copy> Display for Bounds<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use crate::utils::group_digits as g;
        let x = match (self.min, self.max) {
            (None     , None     ) => format!("ALL"),
            (None     , Some(max)) => format!("< {}", g(max)),
            (Some(min), None     ) => format!(">= {}", g(min)),
            (Some(min), Some(max)) => format!("{} <= x < {}", g(min), g(max)),
         };
        f.write_str(&x)
    }
}

impl Display for Iterations {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("number = {}\nsubsets = {}\n", self.number, self.subsets))
    }
}

impl Display for Fov {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("nvoxels = {:?}\nsize = {:?}", self.nvoxels, self.size))
    }
}

impl Display for Scatter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(Bins{bins}) =  self.phi { f.write_str(&format!("phi.bins = {bins}\n"))? };
        if let Some(r)          = &self.r   { f.write_str(&format!(  "phi.r  = {r}\n"   ))? };
        if let Some(z)          = &self.z   { f.write_str(&format!(  "phi.z  = {z}\n"   ))? };
        if let Some(dz)         = &self.dz  { f.write_str(&format!(  "phi.dz = {dz}\n"  ))? };
        if let Some(dt)         = &self.dt  { f.write_str(&format!(  "phi.dt = {dt}\n"  ))? };
        Ok(())
    }
}

impl<T> Display for BinsMax<T>
where
    T: Copy + Debug + FromStr,
    <T as FromStr>::Err: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let BinsMax{ bins, max } = self;
        f.write_fmt(format_args!("{{bins = {bins}, max = {max:?}}}"))
    }
}

impl Display for BinsLength {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let BinsLength{ bins, length } = self;
        f.write_fmt(format_args!("{{bins = {bins}, max = {length:?}}}"))
    }
}
