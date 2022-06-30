//! Configuration file parser for MLEM

use std::fs;
use std::str::FromStr;
use std::path::PathBuf;

use serde::{Deserialize, Deserializer, de};


use crate::{Length, Ratio, Time};

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

fn deserialize_uom_3d_opt<'d, D, T>(deserializer: D) -> Result<Option<(T, T, T)>, D::Error>
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

#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
pub struct Config {

    /// Number of MLEM or OSEM iterations to perform
    #[serde(default = "mandatory")]
    pub iterations: usize,

    /// Number of OSEM subsets per iteration
    #[serde(default = "default_subsets")]
    pub subsets: usize,

    #[serde(default)]
    #[serde(deserialize_with = "deserialize_uom_opt")]
    pub tof: Option<Time>,

    #[serde(default)]
    #[serde(deserialize_with = "deserialize_uom_opt")]
    pub cutoff: Option<Ratio>,

    #[serde(default = "mandatory")]
    pub nvoxels: (usize, usize, usize),

    #[serde(default = "mandatory")]
    #[serde(deserialize_with = "deserialize_uom_3d")]
    pub fov_size: (Length, Length, Length),

    pub scatter: Option<Scatter>,
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

fn default_subsets() -> usize { 1 }

pub fn read_config_file(path: PathBuf) -> Config {
    let config: String = fs::read_to_string(&path)
        .expect(&format!("Couldn't read config file `{:?}`", path));
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

    use geometry::units::{cm, mm, ps};

    // ----- Test an example on-disk config file -----------------------------------------
    #[test]
    fn test_config_file() {
        let config = read_config_file("mlem-config.toml".into());
        assert_eq!(config.iterations, 4);
        assert_eq!(config.subsets, 20);
        assert_eq!(config.tof, Some(ps(200.0)));

        let scatter = config.scatter.unwrap();
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
    #[test]
    fn config_iterations() {
        check!{Config("iterations = 50") fields:
               iterations = 50;
               subsets     = 1
        }

        check!{Config(r#"
                 iterations = 4
                 subsets = 20
               "#) fields:
               iterations =  4;
               subsets    = 20
        }
    }
    // ----- Make sure that unknown fields are not accepted -----------------------------
    #[test]
    #[should_panic]
    fn config_reject_unknown_field() {
        parse("unknown_field = 666")
    }
    // ----- Test FOV parameters ---------------------------------------------------------
    #[test]
    fn config_fov() {
        check!{Config(r#"
                     nvoxels = [10, 20, 30]
                     fov_size = ["123 mm", "456 mm", "78 cm"]
               "#)
        fields:
               nvoxels  = (    10   ,    20    ,    30   );
               fov_size = (mm(123.0), mm(456.0), cm(78.0));
        }
    }
    // ----- Test Scattergram parameters ------------------------------------------------
    #[test]
    fn config_scattergram() {
        let c: Config = parse(r#"
                 [scatter]
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
        let scatter = c.scatter.unwrap();

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
