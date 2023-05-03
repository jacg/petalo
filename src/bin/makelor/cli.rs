/// Command line interface for `makelor` executable
#[derive(clap::Parser, Debug, Clone)]
#[clap(
    name = "makelor",
    about = "Create LORs from MC data",
    subcommand_precedence_over_arg = true, // This can probably be removed in clap 4
)]
pub (super) struct Cli {
    /// HDF5 input files with waveform and charge tables
    pub infiles: Vec<PathBuf>,

    /// HDF5 output file for LORs found in input file
    #[clap(short, long)]
    pub out: PathBuf,

    /// Maximum number of rayon threads used by MLEM
    #[clap(short = 'j', long, default_value = "2")] // HDF reads seems to make > 2 pointless
    pub threads: usize,

    /// Chunk size in output HDF5 file
    #[clap(short = 'c', long, default_value = "1000000")]
    pub chunk_size: usize,

    #[clap(subcommand)]
    pub (super) reco: Reco,

    // TODO allow using different group/dataset in output
}

#[derive(clap::Parser, Debug, Clone)]
pub (super) enum Reco {

    /// Reconstruct LORs from first vertices in LXe
    FirstVertex,

    /// Reconstruct LORs from barycentre of vertices in LXe
    BaryVertex,

    /// Reconstruct LORs from clusters found by splitting cylinder in half
    Half {
        /// Ignore sensors with fewer hits
        #[clap(short, long = "charge-threshold", default_value = "4")]
        q: u32,
    },

    /// Reconstruct LORs from vertices adjusted to element centres
    Discrete {
        /// Inner radius of scintillator
        #[clap(short, long)]
        r_min: Length,

        /// Radial size of elements = thickness of scintillator
        #[clap(long)]
        dr: Length,

        /// Axial size of elements
        #[clap(long)]
        dz: Length,

        /// Azimuthal width of scintillator elements at `r_min + dr/2`?
        #[clap(long)]
        da: Length,
    },

    /// Reconstruct LORs form DBSCAN clusters
    Dbscan {
        /// Ignore sensors with fewer hits
        #[clap(short, long = "charge-threshold", default_value = "4")]
        q: u32,

        /// Minimum number of sensors in cluster
        #[clap(short = 'n', long, default_value = "10")]
        min_count: usize,

        /// Maximum distance between neighbours in cluster
        #[clap(short = 'd', long, default_value = "100 mm")]
        max_distance: Length,
    },

    /// Reconstruct LORs from hits using barycentre of clusters.
    SimpleRec {
        /// Sensor PDE
        #[structopt(short, long, default_value = "0.3")]
        pde: f32,

        /// Sensor time smearing sigma in nanoseconds.
        #[structopt(short, long, default_value = "0.05")]
        sigma_t: f32,

        /// Charge threshold for individual sensors.
        #[structopt(short, long, default_value = "2")]
        threshold: f32,

        /// Minimum number of sensors per cluster.
        #[structopt(short, long, default_value = "2")]
        nsensors: usize,

        /// Charge range to accept sum in clusters.
        #[structopt(short, long, value_parser = parse_bounds::<f32>, default_value = "1..5000")]
        charge_limits: BoundPair<f32>,
    }
}
// ----- Imports -----------------------------------------------------------------------------------------
use std::path::PathBuf;
use units::Length;
use petalo::{
    BoundPair,
    utils::parse_bounds,
};
