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
    #[structopt(long, short="l", default_value = "1000 mm")]
    pub detector_length: Length,

    /// Detector diameter for sensitivity image generation
    #[structopt(long, short="d", default_value = "710 mm")]
    pub detector_diameter: Length,

    /// Number of random LORs to use in sensitivity image generation
    #[structopt(long, short="n", default_value = "5000000")]
    pub n_lors: usize,

    /// Attenuation fiddle factor
    #[structopt(long)]
    pub rho: Lengthf32,

    /// Maximum number of rayon threads
    #[structopt(short = "j", long, default_value = "10")]
    pub num_threads: usize,

}

use structopt::StructOpt;

use std::{error::Error, io::Write};
use std::path::PathBuf;

use petalo::{utils::group_digits, fov::FOV, Lengthf32};
use petalo::image::Image;

use petalo::{Length, Time};
use geometry::uom::ratio;
use petalo::guomc::ConstZero;
use petalo::system_matrix as sm;

fn main() -> Result<(), Box<dyn Error>> {

    let Cli{ input, output, detector_length, detector_diameter, n_lors, rho, num_threads } = Cli::from_args();

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

    let density = Image::from_raw_file(&input)?;
    report_time(&format!("Read density image {:?}", input));


    pre_report(&format!("Creating sensitivity image, using {} LORs ... ", group_digits(n_lors)))?;
    // TODO parallelize
    let lors = find_potential_lors(n_lors, density.fov, detector_length, detector_diameter);
    // To simplify initial Rayon parallelization inside sensitivity_image, pass
    // in a vector, rather than iterator: rayon iter_par over vectors is easier.
    let lors: Vec<_> = lors.collect();

    report_time("\nGenerated LORs");


    let pool = rayon::ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();
    let sensitivity = pool.install(|| Image::sensitivity_image(density.fov, density, &lors, n_lors, rho));
    report_time("done (Calculated sensitivity image)");

    let outfile = output.or_else(|| Some("sensitivity.raw".into())).unwrap();
    std::fs::create_dir_all(PathBuf::from(&outfile).parent().unwrap())?; // TODO turn into utility with cleaner interface
    sensitivity.write_to_raw_file(&outfile)?;
    report_time(&format!("Wrote sensitivity image to {:?}", outfile));
    Ok(())
}

/// Return a vector (size specified in Cli) of LORs with endpoints on cilinder
/// with length and diameter specified in Cli and passing through the FOV
/// specified in Cli.
fn find_potential_lors(n_lors: usize, fov: FOV, detector_length: Length, detector_diameter: Length) -> impl Iterator<Item = sm::LOR> {
    let (l,r) = (detector_length, detector_diameter / 2.0);
    let one_useful_random_lor = move || {
        loop {
            let p1 = random_point_on_cylinder(l, r);
            let p2 = random_point_on_cylinder(l, r);
            if fov.entry(p1, p2).is_some() {
                return Some(petalo::system_matrix::LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0)))
            }
        }
    };
    std::iter::from_fn(one_useful_random_lor).take(n_lors)
}


fn random_point_on_cylinder(l: Length, r: Length) -> petalo::Point {
    use std::f32::consts::TAU;
    use rand::random;
    let z     = l   * (random::<Lengthf32>() - 0.5);
    let theta = TAU *  random::<Lengthf32>();
    let x = r * theta.cos();
    let y = r * theta.sin();
    petalo::Point::new(x, y, z)
}
