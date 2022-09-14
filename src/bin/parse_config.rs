use std::path::PathBuf;
use clap::Parser;

use petalo::config::mlem;

#[derive(clap::Parser, Debug, Clone)]
struct Cli {
    /// Configuration file
    config_file: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    let config = mlem::read_config_file(args.config_file);
    println!("{config:?}");
    Ok(())
}

