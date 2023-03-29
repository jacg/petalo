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
    #[clap(long, short='l')]
    pub detector_length: Length,

    /// Inner radius of scintillator for sensitivity image generation
    #[clap(long, short='r')]
    pub r_min: Length,

    /// Conversion from density to attenuation coefficient in cm^2 / g
    #[clap(long, default_value = "0.095")]
    pub rho_to_mu: Lengthf32,

    /// Maximum number of rayon threads
    #[clap(short = 'j', long, default_value = "30")]
    pub n_threads: usize,

    #[clap(subcommand)]
    detector_type: DetectorType,
}

#[derive(clap::Parser, Debug, Clone)]
enum DetectorType {

    /// Continuous scintillator
    Continuous {
        /// Number of randomly-generated LORs to use
        n_lors: usize,
    },

    /// Discretize scintillator
    Discrete {

        /// Thickness of scintillator tube/ring
        #[clap(long, short = 'r')]
        dr: Length,

        /// Axial width of scintillator elements
        #[clap(long, short = 'z')]
        dz: Length,

        /// Azimuthal width of scintillator elements at `r_min + dr/2`?
        #[clap(long, short = 'a')]
        da: Length,
    },
}

fn main() -> Result<(), Box<dyn Error>> {

    let Cli { input, output, detector_length, r_min, detector_type, rho_to_mu, n_threads } = Cli::parse();

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
    // Convert from [density in kg/m^3] to [mu in mm^-1]
    let attenuation = density_image_into_attenuation_image(density, rho_to_mu);

    // TOF should not be used as LOR attenuation is independent of decay point
    let parameters = Siddon::notof().data();

    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let sensitivity = match detector_type {
        DetectorType::Continuous { n_lors } => {
            pre_report(&format!("Creating sensitivity image, using {} LORs ... ", group_digits(n_lors)))?;
            pool.install(|| continuous::sensitivity_image::<Siddon>(
                detector_length,
                r_min * 2.0,
                parameters,
                &attenuation,
                n_lors))
        },
        DetectorType::Discrete { dr, dz, da } => {
            let discretize = Discretize { dr, dz, da, r_min };
            pool.install(|| discrete::sensitivity_image::<Siddon>(
                detector_length,
                parameters,
                &attenuation,
                discretize,
            ))
        },
    };
    report_time("done");

    let outfile = output.or_else(|| Some("sensitivity.raw".into())).unwrap();
    std::fs::create_dir_all(PathBuf::from(&outfile).parent().unwrap())?; // TODO turn into utility with cleaner interface
    sensitivity.write_to_raw_file(&outfile)?;
    report_time(&format!("Wrote sensitivity image to {:?}", outfile));
    Ok(())
}

/// TODO Just trying an ugly hack for normalizing the image. Do something sensible instead!
fn normalize(data: &mut ImageData, n: usize) { for e in data.iter_mut() { *e /= n as f32 } }

/// Convert from [density in kg/m^3] to [mu in mm^-1]
fn density_image_into_attenuation_image(density: Image, rho_to_mu: AreaPerMass) -> Image {
    let rho_to_mu: f32 = ratio_({
        let kg = kg(1.0);
        let  m = mm(1000.0);
        let rho_unit = kg / (m * m * m);
        let  mu_unit = 1.0 / mm(1.0);
        rho_to_mu / (mu_unit / rho_unit)
    });
    let mut attenuation = density;
    for voxel in &mut attenuation.data {
        *voxel *= rho_to_mu;
    }
    attenuation

}

mod discrete;
mod continuous;

// ----- Imports -----------------------------------------------------------------------------------------
use clap::Parser;

use std::{
    error::Error,
    io::Write,
    path::PathBuf,
};

use petalo::{
    utils::group_digits,
    FOV,
    image::{Image, ImageData},
    system_matrix::{SystemMatrix, Siddon}, discrete::Discretize,
};

use units::{
    Length, AreaPerMass, ratio_, kg, mm, todo::Lengthf32,
};
