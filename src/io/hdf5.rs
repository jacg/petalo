/// Read LORs from HDF5 tables

use std::error::Error;
use std::ops::RangeBounds;
use crate::lorogram::{Scattergram, Prompt};

#[derive(Clone)]
pub struct Args {
    pub input_file: String,
    pub dataset: String,
    pub event_range: Option<std::ops::Range<usize>>,
    pub use_true: bool,
    pub ecut: BoundPair<Energyf32>,
    pub qcut: BoundPair<crate::Chargef32>,
}

use ndarray::{s, Array1};

use crate::{Energyf32, BoundPair};
use crate::Point;
use crate::system_matrix::LOR;

use geometry::uom::{mm, ns, ratio};

pub fn read_table<T: hdf5::H5Type>(filename: &str, dataset: &str, range: Option<std::ops::Range<usize>>) -> hdf5::Result<Array1<T>> {
    let file = ::hdf5::File::open(filename)?;
    let table = file.dataset(dataset)?;
    let data = if let Some(range) = range {
        table.read_slice_1d::<T,_>(s![range])?
    } else {
        table.read_slice_1d::<T,_>(s![..])?
    };
    Ok(data)
}

#[allow(nonstandard_style)]
pub fn read_lors(args: Args, mut scattergram: Option<Scattergram>) -> Result<Vec<LOR>, Box<dyn Error>> {
    let mut rejected = 0;

    // Read LOR data from disk
    let hdf5_lors: Vec<Hdf5Lor> = {
        read_table::<Hdf5Lor>(&args.input_file, &args.dataset, args.event_range.clone())?
            .iter().cloned()
            .filter(|Hdf5Lor{E1, E2, q1, q2, ..}| {
                let eok = args.ecut.contains(E1) && args.ecut.contains(E2);
                let qok = args.qcut.contains(q1) && args.qcut.contains(q2);
                if eok && qok { true }
                else { rejected += 1; false }
            })
            .collect()
    };

    // Fill scattergram, if scatter corrections requested
    if let Some(ref mut scattergram) = &mut scattergram {
        for h5lor @&Hdf5Lor { x1, x2, E1, E2, .. } in &hdf5_lors {
            if x1.is_nan() || x2.is_nan() { continue }
            let prompt = if E1.min(E2) < 511.0 { Prompt::Scatter } else { Prompt::True };
            scattergram.fill(prompt, &LOR::from(h5lor));
        }
    }

    // Convert raw data (Hdf5Lors) to LORs used by MLEM
    let lors: Vec<_> = if let Some(scattergram) = scattergram {
        // With additive correction weights
        hdf5_lors
            .into_iter()
            .map(|hdf5_lor| {
                let mut lor: LOR = hdf5_lor.into();
                lor.additive_correction = scattergram.value(&lor);
                lor
            }).collect()
    } else {
        // Without additive correction weights
        hdf5_lors.into_iter().map(LOR::from).collect()
    };

    let used = lors.len();
    let used_pct = 100 * used / (used + rejected);
    use crate::utils::group_digits as g;
    println!("Using {} LORs (rejected {}, kept {}%)", g(used), g(rejected), used_pct);
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
