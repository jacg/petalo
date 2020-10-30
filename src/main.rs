use std::error::Error;
use structopt::StructOpt;

use petalo::weights::{Point, VoxelBox, Length, LOR};
use petalo::visualize::{lor_weights, Shape};

use petalo::utils::{parse_triplet, parse_lor};

fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();

    let size = args.vbox_size;
    let nvox = args.nvoxels;

    let vbox = VoxelBox::new(size, nvox);
    println!("vbox: {:?}", vbox);

    // LOR on CLI overrides reading LOR from file.
    let lor = if let Some(l) = args.lor {
        l
    } else {
        if let Some(file) = args.file {
            petalo::io::read_lors(file)?[args.event]
        } else {
            *DEFAULT_LOR
        }
    };

    lor_weights(lor, vbox, args.shape, args.threshold, args.tof.map(|x| x / 2.355));

    Ok(())
}

#[macro_use]
extern crate lazy_static;

lazy_static! {
    static ref DEFAULT_LOR: LOR = LOR {
        t1: 10.00,
        t2: 10.03,
        p1: Point::new(242.73334, 0.0, 95.23425),
        p2: Point::new(-236.53258, -54.5143, -94.40064),
    };
}


#[derive(StructOpt, Debug)]
#[structopt(name = "petalo", about = "Visualize LOR interaction with voxels")]
pub struct Cli {

    /// TOF resolution (FWHM) in ps. If not sepcified, TOF is ignored.
    #[structopt(short = "r", long)]
    tof: Option<Length>,

    /// Ignore voxels with weight below this threshold.
    #[structopt(short, long)]
    threshold: Option<Length>,

    /// How to represent voxels. BOX is better for viewing the geometric
    /// weights; BALL is better for viewing TOF weights.
    #[structopt(possible_values = &Shape::variants(), case_insensitive = true, default_value = "box")]
    shape: Shape,

    /// File containing LOR data
    #[structopt(short, long)]
    file: Option<std::path::PathBuf>,

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
    #[structopt(short, long, parse(try_from_str = parse_lor))]
    lor: Option<LOR>,

}
