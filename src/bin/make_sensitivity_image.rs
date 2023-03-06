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
    //let lors = find_potential_lors::continuous(n_lors, density.fov, detector_length, detector_diameter);
    let lors = find_potential_lors::discrete(density.fov, detector_length, detector_diameter);
    report_time(&format!("Generated {} lors", group_digits(lors.len())));
    let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build().unwrap();
    let n_lors = lors.len();
    let job_size = n_lors / n_threads;
    // TOF should not be used as LOR attenuation is independent of decay point
    let parameters = Siddon::notof().data();
    // Convert from [density in kg/m^3] to [mu in mm^-1]
    let attenuation = density_image_into_attenuation_image(density, rho_to_mu);
    let sensitivity = pool.install(|| sensitivity_image::<Siddon, _>(parameters, &attenuation, &lors, job_size));
    report_time("done");

    let outfile = output.or_else(|| Some("sensitivity.raw".into())).unwrap();
    std::fs::create_dir_all(PathBuf::from(&outfile).parent().unwrap())?; // TODO turn into utility with cleaner interface
    sensitivity.write_to_raw_file(&outfile)?;
    report_time(&format!("Wrote sensitivity image to {:?}", outfile));
    Ok(())
}

/// Create sensitivity image by backprojecting `lors` through `attenuation`
/// image.
pub fn sensitivity_image<'l, S, L>(
    parameters : S::Data,
    attenuation: &Image,
    lors       : L,
    job_size   : usize,
) -> Image
where
    S: SystemMatrix,
    L: IntoParallelIterator<Item = &'l LOR>,
    L::Iter: IndexedParallelIterator + Clone,
{
    let lors = lors.into_par_iter();
    let mut backprojection = project_lors::<S,_,_>(parameters, attenuation, lors.clone(), job_size, project_one_lor_sens::<S>);

    // TODO: Just trying an ugly hack for normalizing the image. Do something sensible instead!
    let size = lors.len() as f32;
    for e in backprojection.iter_mut() {
        *e /= size
    }

    Image::new(attenuation.fov, backprojection)
}


mod find_potential_lors {

    use petalo::discrete::Discretize;
    use super::*;

    /// Return a vector (size specified in Cli) of LORs with endpoints on cilinder
    /// with length and diameter specified in Cli and passing through the FOV
    /// specified in Cli.
    pub fn continuous(n_lors: usize, fov: FOV, detector_length: Length, detector_diameter: Length) -> Vec<LOR> {
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
            .collect()
    }

    /// Return a vector of all the LORs that can be made by connecting centres
    /// of the discretized scintillator elements, and which interact with the `fov`.
    pub fn discrete(fov: FOV, detector_length: Length, detector_diameter: Length) -> Vec<LOR> {
        dbg!(detector_length);
        // For prototyping purposes, hard-wire the scintillator element size
        let dz = mm(3.0);
        let da = mm(3.0);
        let dr = mm(30.0);
        let all_centres = Discretize::new(detector_diameter, dr, dz, da)
            .all_element_centres(detector_length)
            .collect::<Vec<_>>();
        dbg!(group_digits(all_centres.len()));

        let origin = petalo::Point::new(mm(0.0), mm(0.0), mm(0.0));
        let all_centres = &all_centres; // Prevent capture by move in closures below
        (0..all_centres.len())
            .into_par_iter()
            .flat_map_iter(|i| (i..all_centres.len())
                           // Rough approximation to 'passes through FOV'
                           //.filter(|&q| origin.distance_to_line(*p, *q) < fov.half_width.z)
                           .filter_map(move |j| {
                               fov.entry(all_centres[i], all_centres[j]).map(|_| {
                                   LOR::new(ns(0.0), ns(0.0), all_centres[i], all_centres[j], ratio(1.0))
                               })
                           })
            )
            .collect()
    }

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
use clap::Parser;

use std::{
    error::Error,
    io::Write,
    path::PathBuf,
};

use petalo::{
    utils::group_digits,
    FOV, LOR,
    image::Image,
    projector::{project_lors, project_one_lor_sens},
    system_matrix::{SystemMatrix, Siddon},
};

use units::{
    Length, Time, AreaPerMass, ratio, ratio_, kg, mm, ns, todo::Lengthf32,
    uom::ConstZero,
};

use rayon::prelude::*;
