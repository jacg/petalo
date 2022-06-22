use std::error::Error;
use std::fs;
use std::path::PathBuf;

use serde::Deserialize;
use structopt::StructOpt;

#[derive(StructOpt, Debug, Clone)]
struct Cli {
    /// Configuration file
    config_file: PathBuf,
}

/// Toy configuration structure for MLEM. Just experimenting for now.
#[derive(Deserialize, Debug)]
pub struct Config {

    /// Number of MLEM or OSEM iterations to perform
    pub iterations: usize,

    /// Number of OSEM subsets per iteration
    pub subsets: usize,

}



fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();


    let config: String = fs::read_to_string(&args.config_file)
        .expect(&format!("Couldn't read config file `{:?}`", args.config_file));

    let config: Config = toml::from_str(&config)?;

    println!("{config:?}");

    Ok(())
}
