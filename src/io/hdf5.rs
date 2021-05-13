/// Read LORs from HDF5 tables

use std::error::Error;

#[derive(Clone)]
pub struct Args {
    pub input_file: String,
    pub dataset: String,
    pub event_range: std::ops::Range<usize>,
    pub use_true: bool,
}

use ndarray::{s, Array1};

use crate::types::{Length, Point};
use crate::weights::{LOR};
type F = Length;

fn read_table<T: hdf5::H5Type>(filename: &str, dataset: &str, range: std::ops::Range<usize>) -> hdf5::Result<Array1<T>> {
    let file = ::hdf5::File::open(filename)?;
    let table = file.dataset(dataset)?;
    let data = table.read_slice_1d::<T,_>(s![range])?;
    Ok(data)
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Event {
    event_id:    f64,
    true_energy: f64,
    true_r1: f64, true_phi1: f64, true_z1: f64, true_t1: f64,
    true_r2: f64, true_phi2: f64, true_z2: f64, true_t2: f64,
    phot_like1: f64, phot_like2: f64,
    reco_r1: f64, reco_phi1: f64, reco_z1: f64, reco_t1: f64,
    reco_r2: f64, reco_phi2: f64, reco_z2: f64, reco_t2: f64,
    not_sel: f64,
}

impl Event {

    fn reco_coords(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let Event{reco_r1, reco_phi1, reco_z1, reco_t1,
                  reco_r2, reco_phi2, reco_z2, reco_t2, ..} = self;
        (*reco_r1, *reco_phi1, *reco_z1, *reco_t1,
         *reco_r2, *reco_phi2, *reco_z2, *reco_t2,)
    }

    fn true_coords(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        let Event{true_r1, true_phi1, true_z1, true_t1,
                  true_r2, true_phi2, true_z2, true_t2, ..} = self;
        (*true_r1, *true_phi1, *true_z1, *true_t1,
         *true_r2, *true_phi2, *true_z2, *true_t2,)
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

pub fn read_lors(args: Args) -> Result<Vec<LOR>, Box<dyn Error>> {
    let it: Vec<LOR> = read_table::<Event>(&args.input_file, &args.dataset, args.event_range.clone())?
        .iter()
        .map(|e| e.to_lor(args.use_true))
        .collect();
    println!("Using {} events", it.len());
    Ok(it)
}


#[cfg(test)]
mod test {

    use super::*;

    #[test] // Test higher-level `read_lors`
    fn read_lors_hdf5() -> hdf5::Result<()> {

        // suppress spamming stdout
        let _suppress_errors = hdf5::silence_errors();

        // First use the reco data to construct the LORs ...
        let args = Args {
            input_file: "src/io/test.h5".into(),
            dataset: "reco_info/table".into(),
            event_range: 0..4,
            use_true: false,
        };
        let lors = read_lors(args.clone()).unwrap();
        assert_eq!(lors[2].p1.coords.x, -120.7552004817734);

        // ... then use the true data.
        let args = Args { use_true: true, ..args };
        let lors = read_lors(args).unwrap();
        assert_eq!(lors[2].p1.coords.x, -120.73839054997248);
        Ok(())
    }

    #[test] // Test lower-level components used by `read_lors`
    fn read_hdf5() -> hdf5::Result<()> {

        let args = Args {
            input_file: "src/io/test.h5".into(),
            dataset: "reco_info/table".into(),
            event_range: 0..4,
            use_true: false,
        };

        let events = read_table::<Event>(&args.input_file, &args.dataset, args.event_range)?;
        assert_eq!(events[2].true_r1, 394.2929992675781);
        assert_eq!(events[2].reco_r1, 394.3750647735979);

        // Leave this here in case we want to regenerate the test file

        // let file = hdf5::File::create("test.h5")?;
        // let reco_info = file.create_group("reco_info")?;
        // let table = reco_info.new_dataset::<Event>().create("table", 10)?;
        // table.write(&events)?;

        Ok(())
    }
}
