// ----------------------------------- CLI -----------------------------------
#[derive(clap::Parser, Debug, Clone)]
#[clap(name = "make_sensitivity_image", about = "Create sensitivity image from a density image")]
pub struct Cli {

    /// The density image of the FOV
    #[clap(short, long)]
    pub input: PathBuf,

    /// Where to write the resulting sensitivity image
    #[clap(short, long)]
    pub output: Option<PathBuf>,

    /// Detector length for sensitivity image generation
    #[clap(long, short='l', default_value = "1000 mm")]
    pub detector_length: Length,

    /// Detector diameter for sensitivity image generation
    #[clap(long, short='d', default_value = "710 mm")]
    pub detector_diameter: Length,

    /// Number of random LORs to use in sensitivity image generation
    #[clap(long, short='n', default_value = "5000000")]
    pub n_lors: usize,

    /// Conversion from density to attenuation coefficient in cm^2 / g
    #[clap(long, default_value = "0.095")]
    pub rho_to_mu: Lengthf32,

    /// Maximum number of rayon threads
    #[clap(short = 'j', long, default_value = "30")]
    pub n_threads: usize,

}

use clap::Parser;

use std::{
    error::Error,
    io::Write,
    path::PathBuf,
};

use petalo::{
    utils::group_digits,
    LOR,
    fov::FOV,
    image::Image,
    mlem::sensitivity_image,
};

use units::{
    Length, Time, AreaPerMass, ratio, kg, mm, todo::Lengthf32,
    uom::ConstZero,
};

fn main() -> Result<(), Box<dyn Error>> {

    let Cli { input, output, detector_length, detector_diameter, n_lors, rho_to_mu, n_threads } = Cli::parse();

    // Interpret rho_to_mu as converting from [rho in g/cm^3] to [mu in cm^-1]
    let rho_to_mu: AreaPerMass = {
        let g = kg(0.001);
        let cm = mm(10.0);
        let rho_unit = g / (cm * cm * cm);
        let  mu_unit = 1.0 / cm;
        rho_to_mu * (mu_unit / rho_unit)
    };

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
    let lors = find_potential_lors(n_lors, density.fov, detector_length, detector_diameter);
    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let sensitivity = pool.install(|| sensitivity_image(density, lors, n_lors, rho_to_mu));
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
fn find_potential_lors(n_lors: usize, fov: FOV, detector_length: Length, detector_diameter: Length) -> impl rayon::iter::ParallelIterator<Item = LOR> {
    let (l,r) = (detector_length, detector_diameter / 2.0);
    let one_useful_random_lor = move |_lor_number| {
        loop {
            let p1 = random_point_on_cylinder(l, r);
            let p2 = random_point_on_cylinder(l, r);
            if fov.entry(p1, p2).is_some() {
                return LOR::new(Time::ZERO, Time::ZERO, p1, p2, ratio(1.0))
            }
        }
    };

    use rayon::prelude::*;
    (0..n_lors)
        .into_par_iter()
        .map(one_useful_random_lor)
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
