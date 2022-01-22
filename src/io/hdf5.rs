/// Read LORs from HDF5 tables

use std::error::Error;
use std::ops::RangeBounds;

#[derive(Clone)]
pub struct Args {
    pub input_file: String,
    pub dataset: String,
    pub event_range: Option<std::ops::Range<usize>>,
    pub use_true: bool,
    pub legacy_input_format: bool,
    pub ecut: BoundPair<Energy>,
    pub qcut: BoundPair<crate::types::Charge>,
}

use ndarray::{s, Array1};

use crate::types::{Length, Point, Energy, BoundPair};
use crate::weights::LOR;
type F = Length;

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

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Event {
    event_id:    f32,
    true_energy: f32,
    true_r1: f32, true_phi1: f32, true_z1: f32, true_t1: f32,
    true_r2: f32, true_phi2: f32, true_z2: f32, true_t2: f32,
    phot_like1: f32, phot_like2: f32,
    reco_r1: f32, reco_phi1: f32, reco_z1: f32, reco_t1: f32,
    reco_r2: f32, reco_phi2: f32, reco_z2: f32, reco_t2: f32,
    not_sel: f32,
}

impl Event {

    fn reco_coords(&self) -> (f32, f32, f32, f32, f32, f32, f32, f32) {
        let &Event{reco_r1, reco_phi1, reco_z1, reco_t1,
                   reco_r2, reco_phi2, reco_z2, reco_t2, ..} = self;
        (reco_r1, reco_phi1, reco_z1, reco_t1,
         reco_r2, reco_phi2, reco_z2, reco_t2,)
    }

    fn true_coords(&self) -> (f32, f32, f32, f32, f32, f32, f32, f32) {
        let &Event{true_r1, true_phi1, true_z1, true_t1,
                  true_r2, true_phi2, true_z2, true_t2, ..} = self;
        (true_r1, true_phi1, true_z1, true_t1,
         true_r2, true_phi2, true_z2, true_t2,)
    }

    fn to_lor(&self, use_true: bool) -> LOR {
        let (r1, phi1, z1, t1,
             r2, phi2, z2, t2) = match use_true {
            true  => self.true_coords(),
            false => self.reco_coords(),
        };
        let x1 = r1 * phi1.cos();
        let y1 = r1 * phi1.sin();

        let x2 = r2 * phi2.cos();
        let y2 = r2 * phi2.sin();

        LOR::new(t1 as F, t2 as F,
                 Point::new(x1 as F, y1 as F, z1 as F),
                 Point::new(x2 as F, y2 as F, z2 as F),
        )
    }

}

#[allow(nonstandard_style)]
pub fn read_lors(args: Args) -> Result<Vec<LOR>, Box<dyn Error>> {
    let mut rejected = 0;
    let it: Vec<LOR> = if ! args.legacy_input_format {
        read_table::<Hdf5Lor>(&args.input_file, &args.dataset, args.event_range.clone())?
            .iter().cloned()
            .filter(|Hdf5Lor{E1, E2, q1, q2, ..}| {
                let eok = args.ecut.contains(&E1) && args.ecut.contains(&E2);
                let qok = args.qcut.contains(&q1) && args.qcut.contains(&q2);
                if eok && qok { true }
                else { rejected += 1; false }
            })
            .map(LOR::from)
            .collect()
    } else {
        read_table::<Event>  (&args.input_file, &args.dataset, args.event_range.clone())?
            .iter()
            .map(|e| e.to_lor(args.use_true))
            .collect()
    };
    let used = it.len();
    let used_pct = 100 * used / (used + rejected);
    use crate::utils::group_digits as g;
    println!("Using {} LORs (rejected {}, kept {}%)", g(used), g(rejected), used_pct);
    Ok(it)
}


#[cfg(test)]
mod test {

    use crate::utils;

    use super::*;
    use assert_approx_eq::assert_approx_eq;

    #[test] // Test higher-level `read_lors`
    fn read_lors_hdf5() -> hdf5::Result<()> {

        // suppress spamming stdout
        let _suppress_errors = hdf5::silence_errors(true);

        // First use the reco data to construct the LORs ...
        let args = Args {
            input_file: "src/io/test.h5".into(),
            dataset: "reco_info/table".into(),
            event_range: Some(0..4),
            use_true: false,
            legacy_input_format: true,
            ecut: utils::parse_bounds("..").unwrap(),
            qcut: utils::parse_bounds("..").unwrap(),
        };
        let lors = read_lors(args.clone()).unwrap();
        assert_approx_eq!(lors[2].p1.coords.x, -120.7552004817734, 1e-5);

        // ... then use the true data.
        let args = Args { use_true: true, ..args };
        let lors = read_lors(args).unwrap();
        assert_approx_eq!(lors[2].p1.coords.x, -120.73839054997248, 1e-5);
        Ok(())
    }

    #[test] // Test lower-level components used by `read_lors`
    fn read_hdf5() -> hdf5::Result<()> {

        let args = Args {
            input_file: "src/io/test.h5".into(),
            dataset: "reco_info/table".into(),
            event_range: Some(0..4),
            use_true: false,
            legacy_input_format: true,
            ecut: utils::parse_bounds("..").unwrap(),
            qcut: utils::parse_bounds("..").unwrap(),
        };

        let events = read_table::<Event>(&args.input_file, &args.dataset, args.event_range)?;
        assert_approx_eq!(events[2].true_r1, 394.2929992675781);
        assert_approx_eq!(events[2].reco_r1, 394.3750647735979);

        // Leave this here in case we want to regenerate the test file

        // hdf5::File::create("test.h5")?
        //     .create_group("reco_info")?
        //     .new_dataset_builder()
        //     .with_data(&events)
        //     .create("table")?;

        Ok(())
    }
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
        Self { dt, p1: Point::new(x1, y1, z1), p2: Point::new(x2, y2, z2)}
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
