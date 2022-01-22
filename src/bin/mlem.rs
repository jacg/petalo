// ----------------------------------- CLI -----------------------------------
use structopt::StructOpt;

use petalo::utils::{parse_triplet, parse_range, parse_bounds, parse_maybe_cutoff, CutoffOption,
                    group_digits};

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "mlem", about = "Maximum Likelyhood Expectation Maximization")]
pub struct Cli {

    /// Number of MLEM iterations to perform
    #[structopt(short, long, default_value = "5")]
    pub iterations: usize,

    /// Field Of View full-widths in mm
    #[structopt(short, long, parse(try_from_str = parse_triplet::<F>), default_value = "300,300,300")]
    pub size: (F, F, F),

    /// Field Of View size in number of voxels
    #[structopt(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "151,151,151")]
    pub nvoxels: (usize, usize, usize),

    /// TOF resolution (sigma) in ps. If not supplied, TOF is ignored
    #[structopt(short = "r", long)]
    pub tof: Option<Time>,

    /// TOF cutoff (✕ sigma). to disable: `-k no` [Rust version only]
    #[structopt(short = "k", default_value = "3", long, parse(try_from_str = parse_maybe_cutoff))]
    pub cutoff: CutoffOption<Ratio>,

    /// Override automatic generation of image output file name
    #[structopt(short, long)]
    pub out_files: Option<String>,

    /// LORs to read in
    #[structopt(short = "f", long, default_value = "data/in/full_body_phantom_reco_combined.h5")]
    pub input_file: String,

    /// The dataset location inside the input file
    #[structopt(short, long, default_value = "reco_info/lors")]
    pub dataset: String,

    /// Which rows of the input file should be loaded
    #[structopt(short, long, parse(try_from_str = parse_range::<usize>))]
    pub event_range: Option<std::ops::Range<usize>>,

    #[cfg(feature = "ccmlem")]
    /// Use the C version of the MLEM algorithm
    #[structopt(short = "c", long)]
    pub use_c: bool,

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
    #[structopt(short = "E", long, parse(try_from_str = parse_bounds::<Energy>), default_value = "..")]
    pub ecut: BoundPair<Energy>,

    /// Ignore events with detected charge/pes outside this range
    #[structopt(short, long, parse(try_from_str = parse_bounds::<Charge>), default_value = "..")]
    pub qcut: BoundPair<Charge>,

}

// --------------------------------------------------------------------------------

use std::error::Error;
use std::path::PathBuf;
use std::fs::create_dir_all;

use petalo::types::{Length, Time, Ratio, Energy, Charge, BoundPair};
use petalo::weights::VoxelBox;
use petalo::mlem::Image;
use petalo::io;

type F = Length;


fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();

    // Set up progress reporting and timing
    use std::time::Instant;
    let mut now = Instant::now();

    let mut report_time = |message: &str| {
        println!("{}: {} ms", message, group_digits(now.elapsed().as_millis()));
        now = Instant::now();
    };

    // Read event data from disk into memory
    let                      Cli{ input_file, dataset, event_range, use_true, legacy_input_format,
                                  ecut, qcut, .. } = args.clone();
    let io_args = io::hdf5::Args{ input_file, dataset, event_range, use_true, legacy_input_format,
                                  ecut, qcut };
    println!("Reading LOR data from disk ...");
    let measured_lors = io::hdf5::read_lors(io_args)?;
    report_time("Loaded LOR data from disk");

    // Define extent and granularity of voxels
    let vbox = VoxelBox::new(args.size, args.nvoxels);
    // TODO: sensitivity matrix, all ones for now
    let sensitivity_matrix = Image::ones(vbox).data;

    let file_pattern = guess_filename(&args);

    // If the directory where results will be written does not exist yet, make it
    create_dir_all(PathBuf::from(format!("{:02}00.raw", file_pattern)).parent().unwrap())?;

    // Perform MLEM iterations
    #[cfg    (feature = "ccmlem") ] let use_c = args.use_c;
    #[cfg(not(feature = "ccmlem"))] let use_c = false;

    if use_c {
        #[cfg(feature = "ccmlem")] run_cmlem(&args, &measured_lors)
    } else {

        #[cfg(not(feature = "serial"))]
        // Set the maximum number of threads used by rayon for parallel iteration
        match rayon::ThreadPoolBuilder::new().num_threads(args.num_threads).build_global() {
            Err(e) => println!("{}", e),
            Ok(_)  => println!("Using up to {} threads.", args.num_threads),
        }

        for (n, image) in (Image::mlem(vbox, &measured_lors, args.tof, args.cutoff, &sensitivity_matrix))
            .take(args.iterations)
            .enumerate() {
                report_time(&format!("Iteration {:2}", n));
                let path = PathBuf::from(format!("{}{:02}.raw", file_pattern, n));
                petalo::io::raw::Image3D::from(&image).write_to_file(&path)?;
                report_time("                               Wrote raw bin");
                // TODO: step_by for print every
            }
    }
    Ok(())
}

fn guess_filename(args: &Cli) -> String {
    if let Some(pattern) = &args.out_files {
        pattern.to_string()
    } else {
        #[cfg    (feature = "ccmlem") ] let c = if args.use_c { "c" } else { "" };
        #[cfg(not(feature = "ccmlem"))] let c = "";
        let (nx, ny, nz) = args.nvoxels;
        let tof = args.tof.map_or(String::from("OFF"), |x| format!("{:.0}", x));
        format!("data/out/{c}mlem/{nx}_{ny}_{nz}_tof_{tof}",
                c=c, nx=nx, ny=ny, nz=nz, tof=tof)
    }
}

// ---- Use the original tofpet3d libmlem (C version), instead of our own Rust version ---

// TODO: this conversion function should really live in the cmlem package, but
// that would require cmlem to depend on the petalo package, because that's
// where types like LOR are defined ... but this crate is currently in the
// petalo package, and it needs to depend on cmlem to call the cmlem function,
// which introduces a circular package dependency, which cargo does not allow.
// The solution is to move this mlem binary crate out of the petalo package, but
// let's just get it working at all, for the time being, and reorganize the
// packages later

#[cfg(feature = "ccmlem")]
use petalo::weights::LOR;

#[cfg(feature = "ccmlem")]
fn run_cmlem(
    args: &Cli,
    lors: &Vec<LOR>
) {
    // Image dimensions
    let (nx, ny, nz) = args.nvoxels;
    let (sx, sy, sz) = args.size;

    // decompose LORs into separate vectors
    let mut x1 = vec![]; let mut y1 = vec![]; let mut z1 = vec![]; let mut t1 = vec![];
    let mut x2 = vec![]; let mut y2 = vec![]; let mut z2 = vec![]; let mut t2 = vec![];
    for lor in lors {
        x1.push(lor.p1.x);
        y1.push(lor.p1.y);
        z1.push(lor.p1.z);
        t1.push(0.0);
        x2.push(lor.p2.x);
        y2.push(lor.p2.y);
        z2.push(lor.p2.z);
        t2.push(petalo::types::ns_to_ps(lor.dt));
    }

    // Add underscore to separate base name from suffix (to match what happens
    // in the Rust version)
    let mut files = guess_filename(&args);
    files.push('_');
    files.push('0'); // Leading zero too!

    // TODO: Dummy sensitivity matrix, for now
    let sensitivity_matrix = vec![1.0; nx * ny * nz];

    cmlem::cmlem(
        args.iterations,
        args.tof.is_some(),
        args.tof.unwrap_or(0.0),
        if sx != sy { panic!("cmlem requires x and y FOVs to be equal") } else { sx },
        sz,
        if nx != ny { panic!("cmlem requires Nx and Ny to be equal") } else { nx },
        nz,
        lors.len(),
        x1, y1, z1, t1,
        x2, y2, z2, t2,
        sensitivity_matrix,
        files,
        1, // save every iteration
    );
}
