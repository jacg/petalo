#[derive(clap::Parser, Debug, Clone)]
#[clap(name = "make_sensitivity_image", about = "Create sensitivity image from a density image")]
pub struct Cli {

    /// The density image to use in the forward projection
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

    /// Output image dimensions, if different from input image (e.g. '10 cm  20 cm  255 mm  30 30 35')
    #[clap(long, value_parser = clap::value_parser!(FOV))]
    pub output_image_dims: Option<FOV>,

    #[clap(subcommand)]
    pub detector_type: DetectorType,
}

#[derive(clap::Parser, Debug, Clone)]
pub enum DetectorType {

    /// Continuous scintillator
    Continuous {
        /// Number of randomly-generated LORs to use
        n_lors: usize,
    },

    /// Discretize scintillator
    Discrete {

        /// Radial size of elements = thickness of scintillator
        #[clap(long)]
        dr: Length,

        /// Axial width of scintillator elements
        #[clap(long)]
        dz: Length,

        /// Azimuthal width of scintillator elements at `r_min + dr/2`?
        #[clap(long)]
        da: Length,

        /// How to adjust the ends of the LORs in the element
        #[clap(long)]
        adjust: Adjust,
    },
}

// ----- Imports -----------------------------------------------------------------------------------------
use std::path::PathBuf;

use petalo::{
    FOV,
    discrete::Adjust,
};
use units::{
    Length, todo::Lengthf32,
};
