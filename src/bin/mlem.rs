// ----------------------------------- CLI -----------------------------------
use structopt::StructOpt;

use petalo::{utils::{parse_triplet, parse_range, parse_bounds, parse_maybe_cutoff, CutoffOption,
                     group_digits}};

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "mlem", about = "Maximum Likelyhood Expectation Maximization")]
pub struct Cli {

    /// Number of MLEM iterations to perform
    #[structopt(short, long, default_value = "5")]
    pub iterations: usize,

    /// Field Of View full-widths in mm
    #[structopt(short, long, parse(try_from_str = parse_triplet::<Length>), default_value = "300mm,300mm,300mm")]
    pub size: (Length, Length, Length),

    /// Field Of View size in number of voxels
    #[structopt(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "151,151,151")]
    pub nvoxels: (usize, usize, usize),

    /// TOF time-resolution sigma (eg '200 ps'). TOF ignored if not supplied
    #[structopt(short, long)]
    pub tof: Option<Time>,

    /// TOF cutoff (✕ sigma). to disable: `-k no` [Rust version only]
    #[structopt(short = "k", default_value = "3", long, parse(try_from_str = parse_maybe_cutoff))]
    pub cutoff: CutoffOption<Ratio>,

    /// Override automatic generation of image output file name
    #[structopt(short, long)]
    pub out_files: Option<String>,

    /// LORs to read in
    #[structopt(short = "f", long, default_value = "MC.h5")]
    pub input_file: String, // TODO replace String with PathBuf here and wherever else appropriate

    /// The dataset location inside the input file
    #[structopt(short, long, default_value = "reco_info/lors")]
    pub dataset: String,

    /// Which rows of the input file should be loaded
    #[structopt(short, long, parse(try_from_str = parse_range::<usize>))]
    pub event_range: Option<std::ops::Range<usize>>,

    /// Sensitivity image to be used for corrections
    #[structopt(long)]
    pub sensitivity_image: Option<PathBuf>,

    /// Use true rather than reco LOR data
    #[structopt(long)]
    use_true: bool,

    #[cfg(not(feature = "serial"))]
    /// Maximum number of rayon threads
    #[structopt(short = "j", long, default_value = "4")]
    pub num_threads: usize,

    /// The input dataset contains the old true/reco r/phi format
    #[structopt(long)]
    pub legacy_input_format: bool,

    /// Ignore events with gamma energy/keV outside this range
    #[structopt(short = "E", long, parse(try_from_str = parse_bounds::<Energyf32>), default_value = "..")]
    pub ecut: BoundPair<Energyf32>,

    /// Ignore events with detected charge/pes outside this range
    #[structopt(short, long, parse(try_from_str = parse_bounds::<Chargef32>), default_value = "..")]
    pub qcut: BoundPair<Chargef32>,

}

// --------------------------------------------------------------------------------

use std::error::Error;
use std::path::PathBuf;
use std::fs::create_dir_all;

use petalo::types::{Energyf32, Chargef32, BoundPair};
use petalo::types::{Length, Time, Ratio};
use petalo::lorogram::Scattergram;
use petalo::weights::{LOR, FOV};
use petalo::mlem::Image;
use petalo::io;
use geometry::uom::{mm_, ratio};


fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();

    // Set up progress reporting and timing
    use std::time::Instant;
    let mut now = Instant::now();

    let mut report_time = |message: &str| {
        println!("{}: {} ms", message, group_digits(now.elapsed().as_millis()));
        now = Instant::now();
    };

    report_time("Startup");

    // Read event data from disk into memory
    let                      Cli{ input_file, dataset, event_range, use_true, legacy_input_format,
                                  ecut, qcut, .. } = args.clone();
    let io_args = io::hdf5::Args{ input_file, dataset, event_range, use_true, legacy_input_format,
                                  ecut, qcut };
    println!("Reading LOR data from disk ...");
    let sgram = io::hdf5::read_scattergram(io_args.clone())?;
    let mut measured_lors = io::hdf5::read_lors(io_args)?;
    report_time("Loaded LOR data from disk");
    calculate_additive_correction(&mut measured_lors, sgram);
    report_time("Added additive correction values to LORs");

    // Define field of view extent and voxelization
    let fov = FOV::new(args.size, args.nvoxels);

    let file_pattern = guess_filename(&args);

    // If the directory where results will be written does not exist yet, make it
    create_dir_all(PathBuf::from(format!("{:02}00.raw", file_pattern)).parent().unwrap())?;

    let sensitivity_image: Option<Image> = args.sensitivity_image.as_ref().map(|path| Image::from_raw_file(path)).transpose()?;
    if let Some(i) = sensitivity_image.as_ref() { assert_image_sizes_match(i, args.nvoxels, args.size) };
    if sensitivity_image.is_some() { report_time("Loaded sensitivity image"); }

    #[cfg(not(feature = "serial"))]
    // Set the maximum number of threads used by rayon for parallel iteration
    match rayon::ThreadPoolBuilder::new().num_threads(args.num_threads).build_global() {
        Err(e) => println!("{}", e),
        Ok(_)  => println!("Using up to {} threads.", args.num_threads),
    }

    for (n, image) in (Image::mlem(fov, &measured_lors, args.tof, args.cutoff, sensitivity_image))
        .take(args.iterations)
        .enumerate() {
            report_time(&format!("Iteration {:2}", n));
            let path = PathBuf::from(format!("{}{:02}.raw", file_pattern, n));
            petalo::io::raw::Image3D::from(&image).write_to_file(&path)?;
            report_time("                               Wrote raw bin");
            // TODO: step_by for print every
        }
    Ok(())
}

fn guess_filename(args: &Cli) -> String {
    if let Some(pattern) = &args.out_files {
        pattern.to_string()
    } else {
        let (nx, ny, nz) = args.nvoxels;
        let tof = args.tof.map_or(String::from("OFF"), |x| format!("{:.0?}", x));
        format!("data/out/mlem/{nx}_{ny}_{nz}_tof_{tof}",
                nx=nx, ny=ny, nz=nz, tof=tof)
    }
}


type FovSize = (Length, Length, Length);
type NVoxels = (usize , usize , usize );

/// Panic if the image size does not match the specified values
fn assert_image_sizes_match(image: &Image, nvoxels: NVoxels, fov_size: FovSize) {
    let size = image.fov.half_width;
    let (idx, idy, idz) = (size[0]*2.0, size[1]*2.0, size[2]*2.0);
    let [inx, iny, inz] = image.fov.n;
    let (enx, eny, enz) = nvoxels;
    let (edx, edy, edz) = fov_size;
    if ! ((enx, eny, enz) == (inx, iny, inz) && (edx, edy, edz) == (idx, idy, idz)) {
        // TODO enable use of density images with different
        // pixelizations as long as they cover the whole FOV.
        println!("Mismatch sensitivity image and output image size:");
        println!("Sensitivity image: {:3} x {:3} x {:3} pixels, {:3} x {:3} x {:3} mm", inx,iny,inz, mm_(idx),mm_(idy),mm_(idz));
        println!("     Output image: {:3} x {:3} x {:3} pixels, {:3} x {:3} x {:3} mm", enx,eny,enz, mm_(edx),mm_(edy),mm_(edz));
        panic!("For now, the sensitivity image must match the dimensions of the output image exactly.");
    }
}

fn calculate_additive_correction(measured_lors: &mut [LOR], sgram: Scattergram) {
    for lor in measured_lors {
        lor.additive_correction = ratio(sgram.value(lor));
    }
}
