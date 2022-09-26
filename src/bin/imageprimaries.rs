use indicatif::{ProgressBar, ProgressStyle};
/// Create raw image from the primary vertices generated in a MC run

// ----------------------------------- CLI -----------------------------------
use clap::Parser;
use std::path::PathBuf;

use petalo::utils::parse_triplet;

#[derive(clap::Parser, Debug, Clone)]
#[clap(setting = clap::AppSettings::ColoredHelp)]
#[clap(name = "imageprimaries", about = "Generate raw image from primary vertices")]
pub struct Cli {

    /// Input file with MC/primaries dataset
    #[clap(short = 'f', long)]
    pub input_files: Vec<PathBuf>,

    /// Image output file
    #[clap(short, long, default_value = "primaries.raw")]
    pub out_file: String,

    /// Field Of View full-widths in mm [default: fit to data]
    #[clap(short, long, parse(try_from_str = parse_triplet::<Length>))]
    pub size: Option<(Length, Length, Length)>,

    /// Field Of View size in number of voxels
    #[clap(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "151,151,151")]
    pub nvoxels: (usize, usize, usize),

}
// --------------------------------------------------------------------------------

use std::error::Error;
use units::{Length, mm, mm_, todo::Lengthf32};
use petalo::fov::FOV;
use petalo::image::Image;
use petalo::io::hdf5::mc::{read_primaries, Primary};
type L = Lengthf32;
use units::uom::si::length::millimeter;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    // --- Process input files -------------------------------------------------------
    let Cli{ input_files, nvoxels, out_file, .. } = args.clone();
    let mut all_events: Vec<Primary> = vec![];
    let mut failed_files = vec![];
    // --- Progress bar --------------------------------------------------------------
    let progress = ProgressBar::new(args.input_files.len() as u64).with_message(args.input_files[0].display().to_string());
    progress.set_style(ProgressStyle::default_bar()
                       .template("Processing file: {msg}\n[{elapsed_precise}] {wide_bar} {pos}/{len} ({eta_precise})")
    );
    progress.tick();
    // --- Collect all events --------------------------------------------------------
    for file in input_files {
        progress.set_message(format!("Processing file {}", file.display()));
        // TODO: replace this loop with a chain of iterators
        if let Ok(events) = read_primaries(&file, petalo::config::mlem::Bounds::none()) {
            for event in events.iter() {
                all_events.push(event.clone());
            }
        } else { failed_files.push(file); }
        progress.inc(1);
    }
    // --- Determine size of output image --------------------------------------------
    let size = if let Some((x,y,z)) = args.size { // from CLI args
        (x,y,z)
    } else { // from extrema of positions in input files
        let (mut xmax, mut ymax, mut zmax): (Length, Length, Length) = (mm(10.0), mm(10.0), mm(10.0));
        for &Primary{ x,y,z, ..} in all_events.iter() {
            xmax = xmax.max(mm(x).abs());
            ymax = ymax.max(mm(y).abs());
            zmax = zmax.max(mm(z).abs());
        }
        ((2.0 * xmax).ceil::<millimeter>(), (2.0 * ymax).ceil::<millimeter>(), (2.0 * zmax).ceil::<millimeter>())
    };
    // --- Create empty image of appropriate size ------------------------------------
    let fov = FOV::new(size, nvoxels);
    let mut image = Image::empty(fov);
    // --- Calculate how to translate spatial position into image index --------------
    let (xn, yn, zn) = args.nvoxels;
    let (xe, ye, ze) = size;
    let (xh, yh, zh) = (xe/2.0, ye/2.0, ze/2.0);
    let (xs, ys, zs) = (mm_(xe) / xn as L,
                        mm_(ye) / yn as L,
                        mm_(ze) / zn as L);
    let pos_to_index3 = |x,y,z| {
        let i = pos_to_index1(x, xs, xn, mm_(xh))?;
        let j = pos_to_index1(y, ys, yn, mm_(yh))?;
        let k = pos_to_index1(z, zs, zn, mm_(zh))?;
        Some([i,j,k])
    };
    // --- Collect event data into image ---------------------------------------------
    for &Primary{ x,y,z, ..} in all_events.iter() {
        if let Some(i3) = pos_to_index3(x,y,z) {
            image[i3] += 1.0;
        }
    }
    // --- Report any files that failed no be read -----------------------------------
    if !failed_files.is_empty() {
        println!("Warning: failed to read the following files:");
        for file in failed_files.iter() {
            println!("  {}", file.display());
        }
        let n = failed_files.len();
        let plural = if n == 1 { "" } else { "s" };
        println!("Warning: failed to read {} file{}:", n, plural);
    }
    // --- Write image to file -------------------------------------------------------

    // TODO: Sometimes the final message fails to appear. Maybe it's being
    // overwritten by the termination of the progress bar. But it also seems
    // that the programs crashes, choking on some file, depending on which other
    // files were processed before it. In some of these cases, the crash report
    // is also hidden.
    let (xe, ye, ze) = (mm_(xe),  mm_(ye),  mm_(ze));
    let message =
        format!("Wrote image with phisical size {} x {} x {} and {} x {} x {} voxels to {}",
                                                xe,  ye,  ze,    xn,  yn,  zn,    out_file);
    progress.finish_with_message(message.clone());
    petalo::io::raw::Image3D::from(&image).write_to_file(out_file)?;
    println!("{}", message);
    Ok(())
}

fn pos_to_index1(position: L, voxel_size: L, nvoxels: usize, half_extent: L) -> Option<usize> {
    if position.abs() >= half_extent { return None }
    let voxel_pos = position / voxel_size;
    let it = ( voxel_pos + (nvoxels as L / 2.0)).floor() as usize;
    Some(it)
}
