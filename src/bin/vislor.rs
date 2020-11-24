use std::error::Error;
use structopt::StructOpt;

use petalo::weights::{VoxelBox, Length, LOR, Ratio};
use petalo::visualize::{lor_weights, Shape};

use petalo::utils::{parse_triplet, parse_lor};
use petalo::io;

fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();

    let size = args.vbox_size;
    let nvox = args.nvoxels;

    let vbox = VoxelBox::new(size, nvox);
    println!("vbox: {:?}", vbox);

    // TODO: reading LOR from file overrides CLI lor: make them mutually
    // exclusive.
    let lor = if let Some(input_file) = args.clone().input_file {
        let event_range = args.event..args.event+1;
        let                      Cli{ dataset, use_true, .. } = args.clone();
        let io_args = io::hdf5::Args{ dataset, use_true, input_file, event_range };
        petalo::io::hdf5::read_lors(io_args)?[0]
    } else {
        args.lor
    };

    lor_weights(lor, vbox, args.shape, args.cutoff, args.sigma);

    Ok(())
}


#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "petalo", about = "Visualize LOR interaction with voxels")]
pub struct Cli {

    /// TOF sensitivity (sigma in ps). If not sepcified, TOF is ignored.
    #[structopt(short = "r", long)]
    sigma: Option<Length>,

    /// Ignore voxels with lie further than <cutoff> from TOF peak.
    #[structopt(short = "k", long)]
    cutoff: Option<Ratio>,

    /// How to represent voxels. BOX is better for viewing the geometric
    /// weights; BALL is better for viewing TOF weights.
    #[structopt(possible_values = &Shape::variants(), case_insensitive = true, default_value = "box")]
    shape: Shape,

    /// LORs to read in
    #[structopt(short = "f", long)]
    pub input_file: Option<String>,

    /// The dataset location inside the input file
    #[structopt(short, long, default_value = "reco_info/table")]
    pub dataset: String,

    /// Event number (in <file>) to be displayed
    #[structopt(short, long, default_value = "0")]
    event: usize,

    /// Dimensions of voxel box in mm
    #[structopt(short, long, parse(try_from_str = parse_triplet::<Length>), default_value = "180,180,180")]
    vbox_size: (Length, Length, Length),

    /// Dimensions of voxel box in voxels
    #[structopt(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "60,60,60")]
    nvoxels: (usize, usize, usize),

    /// LOR to visualize: 't1 t2   x1 y1 z1   x2 y2 z2' (t: ps, xyz: mm)
    #[structopt(short, long, parse(try_from_str = parse_lor), default_value = "0 10  -100 20 -90  100 60 10")]
    lor: LOR,

    /// Use true rather than reco LOR data from file
    #[structopt(long)]
    use_true: bool

}
