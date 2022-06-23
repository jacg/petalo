use std::fs;
use std::path::PathBuf;

use serde::{Deserialize, Deserializer, de};
use structopt::StructOpt;

#[derive(StructOpt, Debug, Clone)]
struct Cli {
    /// Configuration file
    config_file: PathBuf,
}


use petalo::{Length, Time};

fn deserialize_uom<'d, D, T>(deserializer: D) -> Result<Option<T>, D::Error>
where
    D: Deserializer<'d>,
    T: std::str::FromStr,
    <T as std::str::FromStr>::Err: std::fmt::Display,
{
    Option::<&str>::deserialize(deserializer)?
        .map(str::parse::<T>)
        .transpose()
        .map_err(de::Error::custom)
}

/// Toy configuration structure for MLEM. Just experimenting for now.
#[derive(Deserialize, Debug)]
pub struct Config {

    /// Number of MLEM or OSEM iterations to perform
    pub iterations: usize,

    /// Number of OSEM subsets per iteration
    #[serde(default = "default_subsets")]
    pub subsets: usize,

    #[serde(default)]
    #[serde(deserialize_with = "deserialize_uom")]
    pub tof: Option<Time>,

    #[serde(default)]
    #[serde(deserialize_with = "deserialize_uom")]
    pub wtf: Option<Length>,

    //pub nvoxels: (usize, usize, usize),

}

fn default_subsets() -> usize { 1 }

fn read_config_file(path: PathBuf) -> Config {
    let config: String = fs::read_to_string(&path)
        .expect(&format!("Couldn't read config file `{:?}`", path));
    toml::from_str(&config).unwrap()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::from_args();
    let config = read_config_file(args.config_file);
    println!("{config:?}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use geometry::units::{mm, ps};

    // ----- Test an example on-disk config file -----------------------------------------
    #[test]
    fn test_config_file() {
        let config = read_config_file("mlem-config.toml".into());
        assert_eq!(config.iterations, 4);
        assert_eq!(config.subsets, 20);
        assert_eq!(config.tof, Some(ps(200.0)));
        assert_eq!(config.wtf, Some(mm(666.0)));
    }

    // ----- Some helpers to make the tests more concise ---------------------------------
    fn parse<'d, D: Deserialize<'d>>(input: &'d str) -> D {
        toml::from_str(input).unwrap()
    }

    macro_rules! check {
        ($type:ident($text:expr).$field:ident = $expected:expr) => {
            let r: $type = parse::<$type>($text);
            println!("DESERIALIZED: {r:?}");
            assert_eq!(r.$field, $expected);
        };
    }
    // ----- Test deserializing of individual aspects of the Config type ----------------
    #[test]
    fn config_iterations() {
        let mlem = "iterations = 50";
        check!(Config(mlem).iterations = 50);
        check!(Config(mlem).subsets    =  1);

        let osem = r#"
             iterations = 4
             subsets = 20
           "#;
        check!(Config(osem).iterations =  4);
        check!(Config(osem).subsets    = 20);
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
            #[serde(deserialize_with = "deserialize_uom")]
            t: Option<Time>,
            #[serde(default)]
            #[serde(deserialize_with = "deserialize_uom")]
            l: Option<Length>,
        }
        check!(X(r#"t = "2 ps""#).t = Some(ps(2.)));
        check!(X(r#"l = "3 mm""#).l = Some(mm(3.)));
        check!(X(r#"          "#).t = None        );
    }

}
