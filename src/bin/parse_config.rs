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

fn read_config_file(path: PathBuf) -> Config {
    let config: String = fs::read_to_string(&path)
        .expect(&format!("Couldn't read config file `{:?}`", path));
    toml::from_str(&config).unwrap()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let config = read_config_file(args.config_file);
    println!("{config:?}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_file() {
        let config = read_config_file("mlem-config.toml".into());
        assert_eq!(config.iterations, 4);
        assert_eq!(config.subsets, 20);
    }
}
