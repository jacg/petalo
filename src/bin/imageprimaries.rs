/// Create raw image from the primary vertices generated in a MC run

// ----------------------------------- CLI -----------------------------------
use structopt::StructOpt;

use petalo::utils::parse_triplet;

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "imageprimaries", about = "Generate raw image from primary vertices")]
pub struct Cli {

    /// Input file with MC/primaries dataset
    #[structopt(short = "f", long)]
    pub input_files: Vec<String>,

    /// Image output file
    #[structopt(short, long, default_value = "primaries.raw")]
    pub out_file: String,

    /// Field Of View full-widths in mm [default: fit to data]
    #[structopt(short, long, parse(try_from_str = parse_triplet::<L>))]
    pub size: Option<(L, L, L)>,

    /// Field Of View size in number of voxels
    #[structopt(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "151,151,151")]
    pub nvoxels: (usize, usize, usize),

}
// --------------------------------------------------------------------------------

use std::error::Error;
use petalo::types::Length;
use petalo::weights::VoxelBox;
use petalo::mlem::Image;
use petalo::io::hdf5::{read_table, Primary};
type L = Length;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let Cli{ input_files, nvoxels, out_file, .. } = args.clone();
    let mut all_events: Vec<Primary> = vec![];
    let event_range = None;
    let dataset = "MC/primaries";
    for file in input_files {
        // TODO: replace this loop with a chain of iterators
        let events = read_table::<Primary>(&file, dataset, event_range.clone())?; //.iter();
        //let events = &input_files.iter().flat_map(|f| xxx(f));
        for event in events.iter() {
            all_events.push(event.clone());
        }
    }
	  // Determine size of output image
    let size = if let Some((x,y,z)) = args.size { // from CLI args
        (x,y,z)
    } else { // from extrema of positions in input files
        let (mut xmax, mut ymax, mut zmax): (L, L, L) = (10.0, 10.0, 10.0);
        for &Primary{ x,y,z, ..} in all_events.iter() {
            xmax = xmax.max(x.abs());
            ymax = ymax.max(y.abs());
            zmax = zmax.max(z.abs());
        }
        ((2.0 * xmax).ceil(), (2.0 * ymax).ceil(), (2.0 * zmax).ceil())
    };

	// Create empty image of appropriate size
    let vbox = VoxelBox::new(size, nvoxels);
    let mut image = Image::empty(vbox);

    let (xn, yn, zn) = args.nvoxels;
    let (xe, ye, ze) = size;
    let (xh, yh, zh) = (xe/2.0, ye/2.0, ze/2.0);
    let (xs, ys, zs) = (xe / xn as L,
                        ye / yn as L,
                        ze / zn as L);
    let pos_to_index3 = |x,y,z| {
        let i = pos_to_index1(x, xs, xn, xh)?;
        let j = pos_to_index1(y, ys, yn, yh)?;
        let k = pos_to_index1(z, zs, zn, zh)?;
        Some([i,j,k])
    };

    // Collect event data into image
    for &Primary{ x,y,z, ..} in all_events.iter() {
        if let Some(i3) = pos_to_index3(x,y,z) {
            image[i3] += 1.0;
        }
    }

    // Write image to file
    petalo::io::raw::Image3D::from(&image).write_to_file(out_file.clone())?;
    println!("Wrote image with phisical size {} x {} x {} and {} x {} x {} voxels to {}",
                                             xe,  ye,  ze,    xn,  yn,  zn,    out_file);
    Ok(())
}

fn pos_to_index1(position: L, voxel_size: L, nvoxels: usize, half_extent: L) -> Option<usize> {
    if position.abs() >= half_extent { return None }
    let voxel_pos = position / voxel_size;
    let it = ( voxel_pos + (nvoxels as L / 2.0)).floor() as usize;
    Some(it)
}
