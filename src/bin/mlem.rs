// ----------------------------------- CLI -----------------------------------
use structopt::StructOpt;

use petalo::utils::{parse_triplet, parse_range, parse_maybe_cutoff, CutoffOption};

#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "mlem", about = "Maximum Likelyhood Expectation Maximization")]
pub struct Cli {

    /// Number of MLEM iterations to perform
    #[structopt(short, long, default_value = "5")]
    pub iterations: usize,

    /// Voxel box full-widths in mm
    #[structopt(short, long, parse(try_from_str = parse_triplet::<F>), default_value = "180,180,180")]
    pub size: (F, F, F),

    /// Number of voxels in each dimension
    #[structopt(short, long, parse(try_from_str = parse_triplet::<usize>), default_value = "60,60,60")]
    pub n_voxels: (usize, usize, usize),

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
    #[structopt(short, long, default_value = "reco_info/table")]
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

    /// The input dataset contains LORs (not events)
    #[structopt(long)]
    pub read_lors: bool,

}

// --------------------------------------------------------------------------------

use std::error::Error;
use std::path::PathBuf;
use std::fs::create_dir_all;

use petalo::types::{Length, Time, Ratio};
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
        println!("{}: {} ms", message, now.elapsed().as_millis());
        now = Instant::now();
    };

    // Read event data from disk into memory
    let                      Cli{ input_file, dataset, event_range, use_true, read_lors, .. } = args.clone();
    let io_args = io::hdf5::Args{ input_file, dataset, event_range, use_true, read_lors };
    println!("Reading LOR data from disk ...");
    let measured_lors = io::hdf5::read_lors(io_args)?;
    report_time("Loaded LOR data from local disk");

    // Define extent and granularity of voxels
    let vbox = VoxelBox::new(args.size, args.n_voxels);
    // TODO: sensitivity matrix, all ones for now
    let sensitivity_matrix = Image::ones(vbox).data;

    let file_pattern = guess_filename(&args);

    // If the directory where results will be written does not exist yet, make it
    create_dir_all(PathBuf::from(format!("{:02}_00.raw", file_pattern)).parent().unwrap())?;

    // Regions of interest for CRC
    fn polar(r: f32, phi: f32) -> (f32, f32) { (r * phi.cos(), r * phi.sin()) }

    use petalo::fom::ROI;
    let step = std::f32::consts::PI / 6.0;
    let roi_from_centre = 50.0;
    let (hot, cold, bg_activity, bg_radius) = (4.0, 0.0, 1.0, 4.0);
    let rois = vec![
        (ROI::CylinderZ(polar(roi_from_centre,  2.0*step),  4.0),  hot),
        (ROI::CylinderZ(polar(roi_from_centre,  4.0*step),  6.5),  hot),
        (ROI::CylinderZ(polar(roi_from_centre,  6.0*step),  8.5),  hot),
        (ROI::CylinderZ(polar(roi_from_centre,  8.0*step), 11.0),  hot),
        (ROI::CylinderZ(polar(roi_from_centre, 10.0*step), 14.0), cold),
        (ROI::CylinderZ(polar(roi_from_centre, 12.0*step), 18.5), cold),
    ];

    let bg_rois = vec![
        ROI::CylinderZ(polar(roi_from_centre,  1.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  3.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  5.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  7.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  9.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre, 11.0*step), bg_radius),
    ];

    use std::{fs::File, io::{Write, BufWriter}};
    let fom_path = PathBuf::from(format!("{}.crcs", file_pattern));
    let fom_file = File::create(fom_path)?;
    let mut fom_buf = BufWriter::new(fom_file);

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
                let data: Vec<F> = image.data;
                let path = PathBuf::from(format!("{}_{:02}.raw", file_pattern, n));
                write(data.into_iter(), &path)?;
                report_time("  Wrote raw bin");
                // TODO: step_by for print every
                let image = petalo::fom::load_image(&path, vbox)?;
                report_time("  Read raw bin");
                let foms = image.foms(&rois, &bg_rois, bg_activity);
                for crc in &foms.crcs {write!(&mut fom_buf, "{:7.2}", crc)?;}; writeln!(&mut fom_buf)?;
                print!("    CRCs:{:16}",""); for crc in foms.crcs {print!(" {:12.2}", crc);}; println!();
                print!("    SNRs:{:16}",""); for snr in foms.snrs {print!(" {:12.2}", snr);}; println!();
                report_time("  Calculated figures of merit");

            }
    }
    Ok(())
}

fn write(data: impl Iterator<Item = F>, path: &PathBuf) -> Result<(), Box<dyn Error>> {
    use petalo::io::raw::write;
    write(data, path)?;
    Ok(())
}


fn guess_filename(args: &Cli) -> String {
    if let Some(pattern) = &args.out_files {
        pattern.to_string()
    } else {
        #[cfg    (feature = "ccmlem") ] let c = if args.use_c { "c" } else { "" };
        #[cfg(not(feature = "ccmlem"))] let c = "";
        let (nx, ny, nz) = args.n_voxels;
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

use petalo::weights::LOR;

#[cfg(feature = "ccmlem")]
fn run_cmlem(
    args: &Cli,
    lors: &Vec<LOR>
) {
    // Image dimensions
    let (nx, ny, nz) = args.n_voxels;
    let (sx, sy, sz) = args.size;

    // decompose LORs into separate vectors
    let mut x1 = vec![]; let mut y1 = vec![]; let mut z1 = vec![]; let mut t1 = vec![];
    let mut x2 = vec![]; let mut y2 = vec![]; let mut z2 = vec![]; let mut t2 = vec![];
    for lor in lors {
        x1.push(lor.p1.x);
        y1.push(lor.p1.y);
        z1.push(lor.p1.z);
        t1.push(lor.t1);
        x2.push(lor.p2.x);
        y2.push(lor.p2.y);
        z2.push(lor.p2.z);
        t2.push(lor.t2);
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
