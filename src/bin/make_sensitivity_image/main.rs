mod discrete;
mod continuous;
mod cli;
use cli::*;

fn main() -> Result<(), Box<dyn Error>> {

    let Cli { input, output, detector_length, r_min, detector_type, rho_to_mu, n_threads, output_image_dims } = Cli::parse();

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
                n_lors,
                output_image_dims,
            ))
        },
        DetectorType::Discrete { dr, dz, da, adjust } => {
            let discretize = Discretize { dr, dz, da, r_min, adjust };
            pool.install(|| discrete::sensitivity_image::<Siddon>(
                detector_length,
                parameters,
                &attenuation,
                discretize,
                output_image_dims,
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

// ----- Imports -----------------------------------------------------------------------------------------
use std::{
    error::Error,
    io::Write,
    path::PathBuf,
};

use petalo::{
    utils::group_digits,
    FOV,
    image::{Image, ImageData},
    projectors::{SystemMatrix, Siddon}, discrete::Discretize,
};

use units::{
    AreaPerMass, ratio_, kg, mm,
};
use clap::Parser;
