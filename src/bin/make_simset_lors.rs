use std::{path::PathBuf, error::Error, fs::File};
use clap::Parser;
use petalo::io::hdf5::Hdf5Lor;

#[derive(clap::Parser, Debug, Clone)]
#[clap(
    name = "make-simset-lor",
    about = "Create LORs from SimSET history files",
)]
pub struct Args {
    /// SimSET input files with waveform and charge tables
    pub r#in: PathBuf,

    /// HDF5 output file for LORs found in input file
    #[clap(short, long)]
    pub out: PathBuf,

    /// Maximum number of LORs to write
    #[arg(short = 'n', long)]
    stop_after: Option<usize>,

    #[arg(value_enum, short = 't', long = "file-type", default_value_t = HistoryFileType::Custom)]
    history_file_type: HistoryFileType,
}

// TODO: remove copy-paste reuse of enum defined in simset/src/bin/main.rs
#[derive(clap::ValueEnum, Clone, Debug)]
enum HistoryFileType {
    Standard,
    Custom,
}

impl std::fmt::Display for HistoryFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HistoryFileType::Standard => f.write_str("Standard"),
            HistoryFileType::Custom   => f.write_str("Custom"),
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    let mut file = File::open(args.r#in)?;
    simset::skip_header(&mut file)?;
    let lors = std::iter::from_fn(|| Some(lor_from_standard_file(&mut file)));
    for lor in lors.take(args.stop_after.unwrap_or(9999999999999)) {
        println!("{lor:?}");
    }
    Ok(())
}

fn lor_from_standard_file(file: &mut File) -> Hdf5Lor {
    let photons = simset::StandardDecayWithPhotons::read(file).unwrap().photons;
    let mut last_blue = None;
    let mut last_pink = None;
    for current in photons {
        if current.colour() == simset::PhotonColour::Blue { last_blue = Some(current) }
        else                                                 { last_pink = Some(current) }
    }
    let (p1, p2) = (last_blue.unwrap(), last_pink.unwrap());
    Hdf5Lor {
        dt: (p2.time - p1.time) as f32,
        x1: p1.location.x, y1: p1.location.x, z1: p1.location.x,
        x2: p2.location.x, y2: p2.location.y, z2: p2.location.z,
        q1: f32::NAN , q2: f32::NAN,
        E1: p1.energy, E2: p2.energy
    }
}
