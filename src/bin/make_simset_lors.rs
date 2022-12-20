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

    match args.history_file_type {
        HistoryFileType::Standard => {
            simset::standard::iterate_events(&mut file)?
                .take(args.stop_after.unwrap_or(usize::MAX))
                .map(lor_from_standard_event)
                .for_each(show_lor);
        },
        HistoryFileType::Custom => {
            simset::custom::iterate_events(&mut file)?
                .take(args.stop_after.unwrap_or(usize::MAX))
                .filter_map(lor_from_custom_event_maybe)
                .for_each(show_lor);
        },
    }
    Ok(())
}

fn show_lor(Hdf5Lor { dt, x1, y1, z1, x2, y2, z2, E1, E2, .. }: Hdf5Lor) {
    println!("{dt:6.1} ps     ({x1:6.1} {y1:6.1} {z1:6.1})  -mm-  ({x2:6.1} {y2:6.1} {z2:6.1})     {E1:5.1} {E2:5.1}");
}

fn lor_from_standard_event(event: simset::standard::Event) -> Hdf5Lor {
    let mut last_blue = None;
    let mut last_pink = None;
    for current in event.photons {
        if current.colour() == simset::PhotonColour::Blue { last_blue = Some(current) }
        else                                              { last_pink = Some(current) }
    }
    let (p1, p2) = (last_blue.unwrap(), last_pink.unwrap());
    Hdf5Lor {
        dt: s_to_ps(p2.time - p1.time),
        x1: cm_to_mm(p1.x), y1: cm_to_mm(p1.y), z1: cm_to_mm(p1.z),
        x2: cm_to_mm(p2.x), y2: cm_to_mm(p2.y), z2: cm_to_mm(p2.z),
        q1: f32::NAN , q2: f32::NAN,
        E1: p1.energy, E2: p2.energy
    }
}

fn lor_from_custom_event_maybe(simset::custom::Event { blues, pinks }: simset::custom::Event) -> Option<Hdf5Lor> {
    let (p1, p2) = (blues.last()?, pinks.last()?);
    use units::mm_;
    Some(Hdf5Lor {
        dt: cm_to_ps(p2.travel_distance - p1.travel_distance), // TODO check sign
        x1: mm_(p1.x), y1: mm_(p1.y), z1: mm_(p1.z),
        x2: mm_(p2.x), y2: mm_(p2.y), z2: mm_(p2.z),
        q1: f32::NAN, q2: f32::NAN,
        E1: p1.energy as f32, E2: p2.energy as f32,
    })
}

fn  s_to_ps( s: f64) -> f32 {  s as f32 * 10e12 }
fn cm_to_mm(cm: f32) -> f32 { cm * 10.0  }
fn cm_to_ps(cm: f64) -> f32 {
    use units::{C, ps, ratio_};
    let c = ratio_(C / (units::cm(1.0) / ps(1.0)));
    cm as f32 / c
}
