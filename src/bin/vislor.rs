use std::error::Error;
use std::path::PathBuf;
use clap::Parser;

use units::{Time, Ratio, todo::Lengthf32};
use petalo::{system_matrix::LOR, fov::FOV};
use petalo::visualize::{lor_weights, Shape};

use petalo::utils::{parse_triplet, parse_lor, parse_maybe_cutoff, CutoffOption};
use petalo::io;

use petalo::config::mlem::{Config, Bounds, Input};

use units::mm;

fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::parse();

    let (dx, dy, dz) = args.size;
    let size = (mm(dx), mm(dy), mm(dz));
    let nvox = args.nvoxels;

    let fov = FOV::new(size, nvox);
    println!("fov: {fov:?}");

    // TODO: reading LOR from file overrides CLI lor: make them mutually
    // exclusive.
    let lor = if let Some(file) = args.clone().input_file {
        let events = Bounds::<usize> {min: Some(args.event), max: Some(args.event+1)};
        let                      Cli{ dataset, .. } = args.clone();
        let config = Config {
            input: Input { dataset, file, events, ..Default::default()},
            ..Default::default()
        };
        io::hdf5::read_lors(&config, None, 1)?[0]
    } else {
        args.lor
    };

    println!("{}", lor);
    let tof = args.tof.map(|sigma| petalo::config::mlem::Tof { sigma, cutoff: args.cutoff.unwrap() });
    lor_weights(lor, fov, args.shape, tof);
    Ok(())
}


#[derive(Parser, Debug, Clone)]
#[clap(setting = clap::AppSettings::ColoredHelp)]
#[clap(name = "vislor", about = "Visualize LOR interaction with voxels")]
pub struct Cli {

    /// TOF time-resolution sigma (eg '200 ps'). TOF ignored if not supplied
    #[clap(short, long)]
    tof: Option<Time>,

    /// TOF cutoff (✕ sigma). to disable: `-k no`
    #[clap(short = 'k', default_value = "3", long, parse(try_from_str = parse_maybe_cutoff))]
    cutoff: CutoffOption<Ratio>,

    /// How to represent voxels. BOX is better for viewing the geometric
    /// weights; BALL is better for viewing TOF weights.
    #[clap(value_enum, case_insensitive = true, default_value = "box")]
    shape: Shape,

    /// LORs to read in
    #[clap(short = 'f', long)]
    pub input_file: Option<PathBuf>,

    /// The dataset location inside the input file
    #[clap(short, long, default_value = "reco_info/lors")]
    pub dataset: String,

    /// Event number (in <file>) to be displayed
    #[clap(short, long, default_value = "0")]
    event: usize,

    /// Field Of View full-widths in mm
    #[clap(short, long, parse(try_from_str = parse_triplet::<Lengthf32>), default_value = "300,300,300")]
    size: (Lengthf32, Lengthf32, Lengthf32),

    /// Field Of View size in number of voxels
    #[clap(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "151,151,151")]
    nvoxels: (usize, usize, usize),

    /// LOR to visualize: 't1 t2   x1 y1 z1   x2 y2 z2' (t: ps, xyz: mm)
    #[clap(short, long, parse(try_from_str = parse_lor), default_value = "0 300  -100 20 -90  100 60 10")]
    lor: LOR,

}
