/// Read LORs from HDF5 tables

use std::error::Error;
use std::path::Path;
use crate::{
    LOR, Point,
    config::mlem::{Bounds, Config},
    lorogram::{Scattergram, Prompt},
    utils::timing::Progress
};

use rayon::prelude::*;

use ndarray::{s, Array1};

use units::{mm, mm_, ns, ratio, Ratio};


pub fn read_table<T: hdf5::H5Type>(filename: &dyn AsRef<Path>, dataset: &str, events: Bounds<usize>) -> hdf5::Result<Array1<T>> {
    let file = ::hdf5::File::open(filename)?;
    let dataset = file.dataset(dataset)?;
    let Bounds { min, max } = events;
    let data = match (min, max) {
        (None    , None    ) => dataset.read_slice_1d::<T,_>(s![  ..  ])?,
        (Some(lo), None    ) => dataset.read_slice_1d::<T,_>(s![lo..  ])?,
        (None    , Some(hi)) => dataset.read_slice_1d::<T,_>(s![  ..hi])?,
        (Some(lo), Some(hi)) => dataset.read_slice_1d::<T,_>(s![lo..hi])?,
     };
    Ok(data)
}



/// Fill `scattergram`, with spatial distribution of scatters probabilities
/// gathered from `lors`
fn fill_scattergram<'l, L>(
    scattergram: Scattergram,
    lors: L,
    job_size: usize
) -> Scattergram
where
    L: IntoParallelIterator<Item = &'l Hdf5Lor>,
    L::Iter: IndexedParallelIterator,
{
    let empty_scattergram = || scattergram.clone();
    let add_scattergrams = |mut a: Scattergram, b: Scattergram| { a += &b; a };
    let lor_into_scattergram = |mut scattergram: Scattergram, h5lor @&Hdf5Lor { x1, x2, E1, E2, .. }| {
        if x1.is_nan() || x2.is_nan() { return scattergram }
        let prompt = if E1.min(E2) < 510.0 { Prompt::Scatter } else { Prompt::True };
        scattergram.fill(prompt, &LOR::from(h5lor));
        scattergram
    };

    lors.into_par_iter()
        .fold_chunks(job_size, empty_scattergram, lor_into_scattergram)
        .reduce(empty_scattergram, add_scattergrams)

}

/// Read HDF5 LORs from file, potentially filtering according to event, energy
/// and charge ranges
fn read_hdf5_lors(config: &Config) -> Result<(Vec<Hdf5Lor>, usize), Box<dyn Error>> {
    let mut cut = 0;
    let input = &config.input;
    let z_max = config.detector_full_axial_length.map(|l| mm_(l.dz / 2.0));
    // Read LOR data from disk
    let hdf5_lors: Vec<Hdf5Lor> = {
        read_table::<Hdf5Lor>(&input.file, &input.dataset, input.events.clone())?
            .iter().cloned()
            .filter(|&Hdf5Lor{E1, E2, q1, q2, z1, z2, ..}| {
                let eok = input.energy.contains(E1) && input.energy.contains(E2);
                let qok = input.charge.contains(q1) && input.charge.contains(q2);
                let lok = z_max.map_or(true, |z| z1.abs() < z && z2.abs() < z);
                if eok && qok && lok { true }
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
    let (mut hdf5_lors, cut) = read_hdf5_lors(config)?;
    progress.done_with_message(&format!("loaded {}", g(hdf5_lors.len())));

    // Smear energy.
    if let Some(smear_e) = &config.smear_energy { smear_energies(&mut hdf5_lors, smear_e.fwhm, &mut progress); }

    smear_positions(&mut hdf5_lors, config, &mut progress);

    // Use LORs to gather statistics about spatial distribution of scatter probability
    if let Some(sgram) = scattergram {
        progress.start("   Filling scattergram");
        let pool = rayon::ThreadPoolBuilder::new().num_threads(n_threads).build()?;
        let job_size = hdf5_lors.len() / n_threads;
        scattergram = pool.install(|| Some(fill_scattergram(sgram, &hdf5_lors, job_size)));
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

fn smear_energies(hdf5_lors: &mut [Hdf5Lor], fwhm: Ratio, progress: &mut Progress) {
    let smear_energy = |e| {
        use rand_distr::{Normal, Distribution};
        use units::ratio_;
        let fwhm = ratio_(fwhm) * e;
        let sigma = fwhm / 2.35;
        let gauss = Normal::new(e, sigma).unwrap();
        gauss.sample(&mut rand::thread_rng())
    };

    progress.start(&format!("   Smearing energy: {:.0?}% FWHM", fwhm * 100.0));
    for Hdf5Lor { E1, E2, .. } in hdf5_lors.iter_mut() {
        *E1 = smear_energy(*E1);
        *E2 = smear_energy(*E2);
    }
    progress.done();
}

fn smear_positions(hdf5_lors: &mut [Hdf5Lor], config: &Config, progress: &mut Progress) {
    let file = hdf5::File::open(&config.input.file).unwrap();
    let dataset = file.dataset("reco_info/lors").unwrap();
    let get = |attr_name| mm(dataset.attr(attr_name).unwrap()
        .read_1d::<f32>().unwrap()
        .as_slice().unwrap()
        [0]);
    let discretize = crate::discrete::Discretize {
        r_min: get("r_min"),
        dr: get("dr"),
        dz: get("dz"),
        da: get("da"),
        adjust: crate::discrete::Adjust::RandomZPhi,
    };
    let smear_position = discretize.make_adjust_fn();
    progress.start(&format!("   Smearing position: {discretize:.1?}"));
    for Hdf5Lor { x1, y1, z1, x2, y2, z2, .. } in hdf5_lors.iter_mut() {
        let (x,y,z) = smear_position((mm(*x1), mm(*y1), mm(*z1)));
        (*x1, *y1, *z1) = (mm_(x), mm_(y), mm_(z));
        let (x,y,z) = smear_position((mm(*x2), mm(*y2), mm(*z2)));
        (*x2, *y2, *z2) = (mm_(x), mm_(y), mm_(z));
    }
    progress.done();
}

// Include specific table readers and associated types
pub mod mc;
pub mod sensors;



// --------------------------------------------------------------------------------

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

// ----- TESTS ------------------------------------------------------------------------------------------
// Proof of concept: nested compound hdf5 types
#[allow(nonstandard_style)]
#[cfg(test)]
mod test_nested_compound_hdf5 {

    use super::*;

    #[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
    #[repr(C)]
    pub struct Inner {
        pub a: u32,
        pub b: f32,
    }

    #[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
    #[repr(C)]
    pub struct Outer {
        pub id: u32,
        pub r#true: Inner,
        pub   reco: Inner,
    }

    #[test]
    fn roundtrip() -> Result<(), Box<dyn std::error::Error>> {
        let test_data = vec![
            Outer{ id:0, r#true:Inner{ a: 123, b: 4.56 }, reco:Inner{ a: 789, b: 0.12} },
            Outer{ id:1, r#true:Inner{ a: 132, b: 45.6 }, reco:Inner{ a: 798, b: 10.2} },
        ];

        let dir = tempfile::tempdir()?;
        let file_path = dir.path().join("nested-compound.h5");
        let file_path = file_path.to_str().unwrap();

        {
            hdf5::File::create(file_path)?
                .create_group("just-testing")?
                .new_dataset_builder()
                .with_data(&test_data)
                .create("nested")?;
        }
        // read
        let read_data = read_table::<Outer>(&file_path, "just-testing/nested", Bounds::none())?.to_vec();
        assert_eq!(test_data, read_data);
        println!("Test table written to {}", file_path);
        Ok(())
    }
}

// Proof of concept compound HDF5 containing array
#[allow(nonstandard_style)]
#[cfg(test)]
mod test_hdf5_array {

    #[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
    #[repr(C)]
    pub struct Waveform {
        pub event_id: u16,
        pub charge: u8,
        pub ts: [f32; 10],
    }

    impl Waveform {
        fn new(event_id: u16, charge: u8, v: Vec<f32>) -> Self {
            let mut ts: [f32; 10] = [f32::NAN; 10];
            for (from, to) in v.iter().zip(ts.iter_mut()) {
                *to = *from;
            }
            Self { event_id, charge, ts}
        }
    }

    #[test]
    fn testit() -> Result<(), Box<dyn std::error::Error>> {
        let test_data = vec![
            Waveform::new(0, 8, vec![1.2, 3.4]),
            Waveform::new(9, 3, vec![5.6, 7.8, 9.0]),
        ];
        let filename = "just-testing.h5";
        hdf5::File::create(filename)?
            .create_group("foo")?
            .new_dataset_builder()
            .with_data(&test_data)
            .create("bar")?;
        Ok(())
    }
}
