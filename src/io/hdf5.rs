/// Read LORs from HDF5 tables

use std::error::Error;
use std::path::{Path, PathBuf};
use crate::lorogram::{Scattergram, Prompt};
use crate::config::mlem::{Bounds, Config};

use rayon::prelude::*;

#[derive(Clone)]
pub struct Args {
    pub input_file: PathBuf,
    pub dataset: String,
    pub event_range: Option<std::ops::Range<usize>>,
    pub ecut: Bounds<Energyf32>,
    pub qcut: Bounds<Chargef32>,
}

use ndarray::{s, Array1};

use crate::{Point, Energyf32, Chargef32};
use crate::system_matrix::LOR;

use units::{mm, ns, ratio};


pub fn read_table<T: hdf5::H5Type>(filename: &dyn AsRef<Path>, dataset: &str, events: Bounds<usize>) -> hdf5::Result<Array1<T>> {
    let file = ::hdf5::File::open(filename)?;
    let table = file.dataset(dataset)?;
    let Bounds { min, max } = events;
    let data = match (min, max) {
        (None    , None    ) => table.read_slice_1d::<T,_>(s![  ..  ])?,
        (Some(lo), None    ) => table.read_slice_1d::<T,_>(s![lo..  ])?,
        (None    , Some(hi)) => table.read_slice_1d::<T,_>(s![  ..hi])?,
        (Some(lo), Some(hi)) => table.read_slice_1d::<T,_>(s![lo..hi])?,
     };
    Ok(data)
}



/// Fill `scattergram`, with spatial distribution of scatters probabilities
/// gathered from `lors`
fn fill_scattergram(scattergram: Scattergram, lors: &[Hdf5Lor], nthreads: usize) -> Scattergram {
    let empty_scattergram = || scattergram.clone();
    let add_scattergrams = |mut a: Scattergram, b: Scattergram| { a += &b; a };
    let lor_into_scattergram = |mut scattergram: Scattergram, h5lor @&Hdf5Lor { x1, x2, E1, E2, .. }| {
        if x1.is_nan() || x2.is_nan() { return scattergram }
        let prompt = if E1.min(E2) < 510.0 { Prompt::Scatter } else { Prompt::True };
        scattergram.fill(prompt, &LOR::from(h5lor));
        scattergram
    };

    lors.par_iter()
        .with_min_len(lors.len() / nthreads)
        .fold  (empty_scattergram, lor_into_scattergram)
        .reduce(empty_scattergram, add_scattergrams)

}

/// Read HDF5 LORs from file, potentially filtering according to event, energy
/// and charge ranges
fn read_hdf5_lors(config: &Config) -> Result<(Vec<Hdf5Lor>, usize), Box<dyn Error>> {
    let mut cut = 0;
    let input = &config.input;
    // Read LOR data from disk
    let hdf5_lors: Vec<Hdf5Lor> = {
        read_table::<Hdf5Lor>(&input.file, &input.dataset, input.events.clone())?
            .iter().cloned()
            .filter(|&Hdf5Lor{E1, E2, q1, q2, ..}| {
                let eok = input.energy.contains(E1) && input.energy.contains(E2);
                let qok = input.charge.contains(q1) && input.charge.contains(q2);
                if eok && qok { true }
                else { cut += 1; false }
            })
            .collect()
    };
    Ok((hdf5_lors, cut))
}

#[allow(nonstandard_style)]
pub fn read_lors(config: &Config, mut scattergram: Option<Scattergram>, n_threads: usize) -> Result<Vec<LOR>, Box<dyn Error>> {

    let mut progress = crate::utils::timing::Progress::new();

    // Read LORs from file,
    progress.start("   Reading LORs");
    let (hdf5_lors, cut) = read_hdf5_lors(config)?;
    progress.done_with_message(&format!("loaded {}", g(hdf5_lors.len())));

    // Use LORs to gather statistics about spatial distribution of scatter probability
    if let Some(sgram) = scattergram {
        progress.start("   Filling scattergram");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build()?;
        scattergram = pool.install(|| Some(fill_scattergram(sgram, &hdf5_lors, n_threads)));
        progress.done();
    }
    progress.start("   Converting HDF5 LORs into MLEM LORs");

    let hdf5lor_to_lor: Box<dyn Fn(Hdf5Lor) -> LOR + Sync> = if let Some(scattergram) = scattergram.as_ref() {
        Box::new(|hdf5_lor: Hdf5Lor| {
            let mut lor: LOR = hdf5_lor.into();
            lor.additive_correction = scattergram.value(&lor);
            lor
        })
    } else { Box::new(LOR::from) };

    // Convert raw data (Hdf5Lors) to LORs used by MLEM
    let lors: Vec<_> = hdf5_lors
        // This seems to scale pretty perfectly with number of threads, so we
        // let Rayon dynamically pick the actual number of threads, so we don't
        // run this in a local thread pool.
        .into_par_iter()
        .map(&hdf5lor_to_lor)
        .collect();
    progress.done();

    let used = lors.len();
    let used_pct = 100 * used / (used + cut);
    use crate::utils::group_digits as g;
    println!("   Using {} LORs (cut {}    kept {}%)",
                  g(used),      g(cut),   used_pct);
    Ok(lors)
}



// --------------------------------------------------------------------------------
#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct SensorXYZ {
    pub sensor_id: u32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Charge {
    pub event_id: u64,
    pub sensor_id: u64,
    pub charge: u64,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct SensorHit {
    pub event_id: u64,
    pub sensor_id: u64,
    pub time: f64,
}

// The LOR used by mlem contains fields (the points) with types (ncollide Point)
// which hdf5 appears not to be able to digest, so hack around the problem for
// now, by creating a LOR type that is hdf5able.
// Additionally, we now want to store extra information corresponding to the LOR
// (energies, charges) which are useful for applying different cuts later on,
// but irrelevant to MLEM, so two separate LOR types (with and without metadata)
// might actually be the right way to go.
#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
#[allow(nonstandard_style)]
pub struct Hdf5Lor {
    pub dt: f32,
    pub x1: f32,
    pub y1: f32,
    pub z1: f32,
    pub x2: f32,
    pub y2: f32,
    pub z2: f32,
    pub q1: f32,
    pub q2: f32,
    pub E1: f32,
    pub E2: f32,
}

impl From<Hdf5Lor> for LOR {
    fn from(lor: Hdf5Lor) -> Self {
        let Hdf5Lor{dt, x1, y1, z1, x2, y2, z2, ..} = lor;
        Self {
            dt: ns(dt),
            p1: Point::new(mm(x1), mm(y1), mm(z1)),
            p2: Point::new(mm(x2), mm(y2), mm(z2)),
            additive_correction: ratio(1.0)
        }
    }
}

impl From<&Hdf5Lor> for LOR {
    fn from(lor: &Hdf5Lor) -> Self {
        let &Hdf5Lor{dt, x1, y1, z1, x2, y2, z2, ..} = lor;
        Self {
            dt: ns(dt),
            p1: Point::new(mm(x1), mm(y1), mm(z1)),
            p2: Point::new(mm(x2), mm(y2), mm(z2)),
            additive_correction: ratio(1.0)
        }
    }
}

// --------------------------------------------------------------------------------
#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Primary {
    pub event_id: u32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub vx: f32,
    pub vy: f32,
    pub vz: f32,
}
