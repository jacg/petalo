use std::path::PathBuf;
use structopt::StructOpt;

use petalo::config::mlem;

#[derive(StructOpt, Debug, Clone)]
struct Cli {
    /// Configuration file
    config_file: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::from_args();
    let config = mlem::read_config_file(args.config_file);
    println!("{config:?}");
    Ok(())
}

