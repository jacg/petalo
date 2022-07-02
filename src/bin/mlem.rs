use petalo::config::mlem::AttenuationCorrection as AC;
// ----------------------------------- CLI -----------------------------------
use structopt::StructOpt;

use petalo::{lorogram::BuildScattergram, config};

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "mlem", about = "Maximum Likelyhood Expectation Maximization")]
pub struct Cli {

    /// MLEM config file
    pub config_file: PathBuf,

    /// Directory in which results should be written
    pub output_directory: String,

    /// Maximum number of rayon threads
    #[structopt(short = "j", long, default_value = "4")]
    pub num_threads: usize,

}

// --------------------------------------------------------------------------------

use std::error::Error;
use std::path::PathBuf;
use std::fs::create_dir_all;

use petalo::Length;
use petalo::lorogram::Scattergram;
use petalo::fov::FOV;
use petalo::image::Image;
use petalo::io;
use petalo::utils::timing::Progress;
use geometry::units::mm_;


fn main() -> Result<(), Box<dyn Error>> {

    let args = Cli::from_args();
    let config = config::mlem::read_config_file(args.config_file.clone());

    // Set up progress reporting and timing
    let mut progress = Progress::new();

    // Check that output directory is writable. Do this *before* expensive
    // setup, so it fails early
    // If the directory where results will be written does not exist yet, make it
    create_dir_all(PathBuf::from(format!("{:02}00.raw", args.output_directory)).parent().unwrap())
        .expect(&format!("Cannot write in output directory `{}`", args.output_directory));

    // Define field of view extent and voxelization
    let fov = FOV::new(config.fov.size, config.fov.nvoxels);

    let scattergram = build_scattergram(&config.scatter_correction);
    progress.done_with_message("Startup");

    let sensitivity_image =
        if let Some(AC { sensitivity_image: path }) = config.attenuation_correction.as_ref() {
            let image = Image::from_raw_file(&path)
                .expect(&format!("Cannot read sensitivity image {:?}", path.display()));
            assert_image_sizes_match(&image, config.fov.nvoxels, config.fov.size);
            progress.done_with_message("Loaded sensitivity image");
            Some(image)
        } else { None };

    progress.startln("Loading LORs from file");
    let measured_lors = io::hdf5::read_lors(&config, scattergram)?;
    progress.done_with_message("Loaded LORs from file");

    // Set the maximum number of threads used by rayon for parallel iteration
    match rayon::ThreadPoolBuilder::new().num_threads(args.num_threads).build_global() {
        Err(e) => println!("{}", e),
        Ok(_)  => println!("Using up to {} threads.", args.num_threads),
    }

    for (image, iteration, subset) in (Image::mlem(fov, &measured_lors, config.tof, sensitivity_image, config.iterations.subsets))
        .take(config.iterations.number * config.iterations.subsets) {
            progress.done_with_message(&format!("Iteration {iteration:2}-{subset:02}"));
            let path = PathBuf::from(format!("{}{iteration:02}-{subset:02}.raw", args.output_directory));
            petalo::io::raw::Image3D::from(&image).write_to_file(&path)?;
            progress.done_with_message("                               Wrote raw bin");
            // TODO: step_by for print every
        }
    Ok(())
}


type FovSize = (Length, Length, Length);
type NVoxels = (usize , usize , usize );

/// Panic if the image size does not match the specified values
fn assert_image_sizes_match(image: &Image, nvoxels: NVoxels, fov_size: FovSize) {
    use float_eq::float_eq;
    let size = image.fov.half_width;
    let (idx, idy, idz) = (size[0]*2.0, size[1]*2.0, size[2]*2.0);
    let [inx, iny, inz] = image.fov.n;
    let (enx, eny, enz) = nvoxels;
    let (edx, edy, edz) = fov_size;
    // Unwrap uom, to make float_eq! work
    let ids = [mm_(idx), mm_(idy), mm_(idz)];
    let eds = [mm_(edx), mm_(edy), mm_(edz)];
    if ! ((enx, eny, enz) == (inx, iny, inz) && float_eq!(eds, ids, ulps_all <= 1)) {
        // TODO enable use of density images with different
        // pixelizations as long as they cover the whole FOV.
        println!("Mismatch sensitivity image and output image size:");
        println!("Sensitivity image: {:3} x {:3} x {:3} pixels, {:3} x {:3} x {:3} mm", inx,iny,inz, mm_(idx),mm_(idy),mm_(idz));
        println!("     Output image: {:3} x {:3} x {:3} pixels, {:3} x {:3} x {:3} mm", enx,eny,enz, mm_(edx),mm_(edy),mm_(edz));
        panic!("For now, the sensitivity image must match the dimensions of the output image exactly.");
    }
}

// TODO: Make builder use config::mlem::Scatter directly
fn build_scattergram(scatter: &Option<config::mlem::Scatter>) -> Option<Scattergram> {
    use petalo::config::mlem::{Bins, BinsMax, BinsLength};
    if let Some(scatter) = scatter.as_ref() {
        let mut builder = BuildScattergram::new();

        if let Some(Bins { bins }) = scatter.phi {
            builder = builder.phi_bins(bins);
        }
        if let Some(BinsMax { bins, max }) = scatter.r {
            builder = builder.r_bins(bins);
            builder = builder.r_max (max );
        }
        if let Some(BinsMax { bins, max }) = scatter.dz {
            builder = builder.dz_bins(bins);
            builder = builder.dz_max (max );
        }
        if let Some(BinsMax { bins, max }) = scatter.dt {
            builder = builder.dt_bins(bins);
            builder = builder.dt_max (max );
        }
        if let Some(BinsLength { bins, length }) = scatter.z {
            builder = builder.z_bins  (bins);
            builder = builder.z_length(length);
        }
        builder.build()
    } else {
        None
    }
}
