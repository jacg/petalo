use std::error::Error;
use structopt::StructOpt;

use petalo::types::{Length, Ratio};
use petalo::weights::{VoxelBox, LOR};
use petalo::visualize::{lor_weights, Shape};

use petalo::utils::{parse_triplet, parse_lor, parse_maybe_cutoff, parse_bounds, CutoffOption};
use petalo::io;

fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();

    let size = args.size;
    let nvox = args.nvoxels;

    let vbox = VoxelBox::new(size, nvox);
    println!("vbox: {:?}", vbox);

    // TODO: reading LOR from file overrides CLI lor: make them mutually
    // exclusive.
    let lor = if let Some(input_file) = args.clone().input_file {
        let event_range = args.event..args.event+1;
        let                      Cli{ dataset, use_true, legacy_input_format, .. } = args.clone();
        let io_args = io::hdf5::Args{ dataset, use_true, legacy_input_format, input_file,
                                      ecut: parse_bounds("..").unwrap(), qcut: parse_bounds("..").unwrap(),
                                      event_range: Some(event_range) };
        petalo::io::hdf5::read_lors(io_args)?[0]
    } else {
        args.lor
    };

    println!("{}", lor);
    lor_weights(lor, vbox, args.shape, args.cutoff, args.sigma);
    Ok(())
}


#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "vislor", about = "Visualize LOR interaction with voxels")]
pub struct Cli {

    /// TOF sensitivity (sigma in ps). If not sepcified, TOF is ignored.
    #[structopt(short = "r", long)]
    sigma: Option<Length>,

    /// TOF cutoff (✕ sigma). to disable: `-k no`
    #[structopt(short = "k", default_value = "3", long, parse(try_from_str = parse_maybe_cutoff))]
    cutoff: CutoffOption<Ratio>,

    /// How to represent voxels. BOX is better for viewing the geometric
    /// weights; BALL is better for viewing TOF weights.
    #[structopt(possible_values = &Shape::variants(), case_insensitive = true, default_value = "box")]
    shape: Shape,

    /// LORs to read in
    #[structopt(short = "f", long)]
    pub input_file: Option<String>,

    /// The dataset location inside the input file
    #[structopt(short, long, default_value = "reco_info/lors")]
    pub dataset: String,

    /// Event number (in <file>) to be displayed
    #[structopt(short, long, default_value = "0")]
    event: usize,

    /// Field Of View full-widths in mm
    #[structopt(short, long, parse(try_from_str = parse_triplet::<Length>), default_value = "300,300,300")]
    size: (Length, Length, Length),

    /// Field Of View size in number of voxels
    #[structopt(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "151,151,151")]
    nvoxels: (usize, usize, usize),

    /// LOR to visualize: 't1 t2   x1 y1 z1   x2 y2 z2' (t: ps, xyz: mm)
    #[structopt(short, long, parse(try_from_str = parse_lor), default_value = "0 300  -100 20 -90  100 60 10")]
    lor: LOR,

    /// Use true rather than reco LOR data from file
    #[structopt(long)]
    use_true: bool,

    /// The input dataset contains the old true/reco r/phi format
    #[structopt(long)]
    legacy_input_format: bool,

}
