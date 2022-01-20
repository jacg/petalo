// ----------------------------------- CLI -----------------------------------
#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "make_sensitivity_image", about = "Create sensitivity image from a density image")]
pub struct Cli {

    /// The density image of the FOV
    #[structopt(short, long)]
    pub input: PathBuf,

    /// Where to write the resulting sensitivity image
    #[structopt(short, long)]
    pub output: Option<PathBuf>,

    /// Detector length for sensitivity image generation
    #[structopt(long, short="l", default_value = "1000")]
    pub detector_length: Length,

    /// Detector diameter for sensitivity image generation
    #[structopt(long, short="d", default_value = "710")]
    pub detector_diameter: Length,

    /// Number of random LORs to use in sensitivity image generation
    #[structopt(long, short="n", default_value = "5000000")]
    pub n_lors: usize,

    /// Attenuation fiddle factor
    #[structopt(long)]
    pub rho: Length,

}

use structopt::StructOpt;

use std::{error::Error, io::Write};
use std::path::PathBuf;

use petalo::{mlem, utils::group_digits, weights::VoxelBox, types::Length};

fn main() -> Result<(), Box<dyn Error>> {

    let Cli{ input, output, detector_length, detector_diameter, n_lors, rho } = Cli::from_args();

    // Set up progress reporting and timing
    use std::time::Instant;
    let mut now = Instant::now();

    let mut report_time = |message: &str| {
        println!("{}: {} ms", message, group_digits(now.elapsed().as_millis()));
        now = Instant::now();
    };

    let pre_report = |message: &str| {
        print!("{}", message);
        std::io::stdout().flush()
    };

    report_time("Startup");

    let density = mlem::Image::from_raw_file(&input)?;
    report_time(&format!("Read density image {:?}", input));

    pre_report(&format!("Generating {} LORs for sensitivity image projections ... ", group_digits(n_lors)))?;
    let lors = find_potential_lors(n_lors, density.vbox, detector_length, detector_diameter);
    report_time("done");

    pre_report("Creating sensitivity image ... ")?;
    // TODO parallelize
    // TODO try to do it without holding all LORs in memory at once (good image uses > 25G RAM!)
    let sensitivity = mlem::Image::sensitivity_image(density.vbox, density, &lors, rho);
    report_time("done");

    let outfile = output.or_else(|| Some("sensitivity.raw".into())).unwrap();
    std::fs::create_dir_all(PathBuf::from(&outfile).parent().unwrap())?; // TODO turn into utility with cleaner interface
    sensitivity.write_to_raw_file(&outfile)?;
    report_time(&format!("Wrote sensitivity image to {:?}", outfile));
    Ok(())
}

/// Return a vector (size specified in Cli) of LORs with endpoints on cilinder
/// with length and diameter specified in Cli and passing through the FOV
/// specified in Cli.
fn find_potential_lors(n_lors: usize, fov: VoxelBox, detector_length: Length, detector_diameter: Length) -> Vec<petalo::weights::LOR> {
    let mut lors = Vec::with_capacity(n_lors);
    let (l,r) = (detector_length, detector_diameter / 2.0);
    while lors.len() < n_lors {
        let p1 = random_point_on_cylinder(l,r);
        let p2 = random_point_on_cylinder(l,r);
        if let Some(_) = fov.entry(&p1, &p2) {
            lors.push(petalo::weights::LOR::new(0.0, 0.0, p1, p2))
        }
    }
    lors
}


fn random_point_on_cylinder(l: Length, r: Length) -> petalo::types::Point {
    use std::f32::consts::TAU;
    use rand::random;
    let z     = l   * (random::<Length>() - 0.5);
    let theta = TAU *  random::<Length>();
    let x = r * theta.cos();
    let y = r * theta.sin();
    petalo::types::Point::new(x, y, z)
}
