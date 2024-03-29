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

use units::{mm, mm_, ns, ratio, ratio_, Ratio};


pub fn read_dataset<T: hdf5::H5Type>(filename: &dyn AsRef<Path>, dataset: &str, events: Bounds<usize>) -> hdf5::Result<Array1<T>> {
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

// Plan
// 1. DONE Assume (panic) chunked dataset, iterate over all items in all chunks
// 2. DONE Add bounds
// 3. TODO Cater for unchunked version
pub fn iter_dataset<T: hdf5::H5Type>(
    filename: &dyn AsRef<Path>,
    dataset: &str,
    Bounds { min, max }: Bounds<usize>,
) -> hdf5::Result<Box<dyn Iterator<Item = T>>>
{
    let file = ::hdf5::File::open(filename)?;
    let dataset = file.dataset(dataset)?;

    let dataset_shape = dataset.shape();
    assert_eq!(dataset_shape.len(), 1);              // Assuming 1-D dataset
    let chunk_dimensions = dataset.chunk().unwrap(); // Assuming dataset is chunked, for now
    assert_eq!(chunk_dimensions.len(), 1);           // Assuming 1-D chunks
    let chunk_size = chunk_dimensions[0];
    let dataset_size = dataset_shape[0];

    let start = min.map_or(0,            |lo| (lo / chunk_size    ) * chunk_size);
    let skip  = min.map_or(0,            |lo| (lo % chunk_size    )             );
    let stop  = max.map_or(dataset_size, |hi| (hi / chunk_size + 1) * chunk_size);
    let take  = max.unwrap_or(dataset_size) - min.unwrap_or(0);

    Ok(Box::new(
        (start..stop)
            .step_by(chunk_size)
            .map(move |n| (n, n+chunk_size))
            .flat_map(move |(b,e)| {
                dataset.read_slice_1d::<T, _>(s![b..e]).into_iter().flatten()
            })
            .skip(skip)
            .take(take)
    ))
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
fn read_hdf5_lors(config: &Config) -> Result<Vec<Hdf5Lor>, Box<dyn Error>> {
    let input = &config.input;
    let z_max = config.detector_full_axial_length.map(|l| mm_(l.dz / 2.0));
    let total = ::hdf5::File::open(&config.input.file).unwrap()
        .dataset(&config.input.dataset).unwrap()
        .shape()[0];
    let Bounds { min, max } = config.input.events;
    let to_be_read = max.map_or(total, |max| max.min(total)) -
                     min.map_or(0    , |min| min.max(0    ));

    let smear_energy = make_smear_energy(config.smear_energy.unwrap().fwhm);
    let pre_smearing_e_cut = 511.0 * (1.0 - 1.6 * ratio_(config.smear_energy.unwrap().fwhm));
    println!("Applying pre-cut at {pre_smearing_e_cut} keV, final cut at {:?}", config.input.energy);

    let progress = progress::Progress::new(to_be_read);
    // Read LOR data from disk
    let hdf5_lors: Vec<Hdf5Lor> = {
        iter_dataset::<Hdf5Lor>(&input.file, &input.dataset, input.events.clone())?
            .inspect(|_|  { progress.read() })
            .filter(|&Hdf5Lor{z1, z2, ..}| { z_max.map_or(true, |z| z1.abs() < z && z2.abs() < z) })
            .filter(|&Hdf5Lor{q1, q2, ..}| { input.charge.contains(q1) && input.charge.contains(q2) })
            .filter(|&Hdf5Lor{E1, E2, ..}| {  pre_smearing_e_cut < E1  &&  pre_smearing_e_cut < E2 })
            .map(|mut l@Hdf5Lor { E1, E2, .. }| {
                l.E1 = smear_energy(E1);
                l.E2 = smear_energy(E2);
                l
            })
            .filter(|&Hdf5Lor{E1, E2, ..}| { input.energy.contains(E1) && input.energy.contains(E2) })
            .inspect(|_|  { progress.selected() })
            .collect()
    };
    progress.done();
    Ok(hdf5_lors)
}

mod progress {
    use std::sync::Mutex;
    use indicatif::{ProgressBar, ProgressStyle};

    use crate::utils::group_digits;

    pub (super) struct Progress(Mutex<Inner>);

    struct Inner {
        read: u64,
        selected: u64,
        bar: ProgressBar,
    }

    impl Progress {

        pub (super) fn new(total: usize) -> Self {
            let bar = ProgressBar::new(total as u64).with_message("Reading and filtering LORs");
            bar.set_style(
                ProgressStyle::default_bar()
                    .template("{msg}\n[{elapsed_precise}] {wide_bar} {human_pos} / {human_len} ({eta_precise})")
                    .unwrap()
            );
            bar.tick();
            Self(Mutex::new(Inner { read: 0, selected: 0, bar }))
        }

        pub (super) fn read(&self) {
            let mut i = self.0.lock().unwrap(); i.read     += 1;
            if i.read % 100000 == 0 {
                i.bar.set_position(i.read);
                     let percent = (100 * i.selected) as f32 / i.read as f32;
                i.bar.set_message(format!(
                    "Kept {} out of {} LORs ({percent:.1}%)", group_digits(i.selected), group_digits(i.read),
                ));
            }
        }
        pub (super) fn selected(&self) {
            let mut i = self.0.lock().unwrap();
            i.selected += 1;
        }
        pub (super) fn done(&self) {
            let i = self.0.lock().unwrap();
            let percent = (100 * i.selected) as f32 / i.read as f32;
            i.bar.finish_with_message(format!(
                "Kept {} out of {} LORs ({percent:.1}%)", group_digits(i.selected), group_digits(i.read),
            ));
            //i.bar.finish();
        }
    }
}


#[allow(nonstandard_style)]
pub fn read_lors(config: &Config, mut scattergram: Option<Scattergram>, n_threads: usize) -> Result<Vec<LOR>, Box<dyn Error>> {

    let mut progress = crate::utils::timing::Progress::new();

    // Read LORs from file,
    progress.start("   Reading LORs");
    let mut hdf5_lors = read_hdf5_lors(config)?;
    use crate::utils::group_digits as g;
    progress.done_with_message(&format!("loaded {}", g(hdf5_lors.len())));

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


    Ok(lors)
}

fn make_smear_energy(fwhm: Ratio) -> impl Fn(f32) -> f32 {
    move |e| {
        use rand_distr::{Normal, Distribution};
        let fwhm = ratio_(fwhm) * e;
        let sigma = fwhm / 2.35;
        let gauss = Normal::new(e, sigma).unwrap();
        gauss.sample(&mut rand::thread_rng())
    }
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

    hdf5_lors.par_iter_mut().for_each(|Hdf5Lor { x1, y1, z1, x2, y2, z2, .. }| {
        let (x,y,z) = smear_position((mm(*x1), mm(*y1), mm(*z1)));
        (*x1, *y1, *z1) = (mm_(x), mm_(y), mm_(z));
        let (x,y,z) = smear_position((mm(*x2), mm(*y2), mm(*z2)));
        (*x2, *y2, *z2) = (mm_(x), mm_(y), mm_(z));
    });

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
#[cfg(test)]
mod test_chunked {
    use super::*;
    use rstest::rstest;
    use pretty_assertions::assert_eq;
    use tempfile::tempdir;
    use itertools::Itertools;

    #[rstest(
        case::front_back_none_back(12, 4, None   , None    , 0..12),
        case::front_part_back_none(12, 4, Some(1), None    , 1..12),
        case::front_full_back_none(12, 4, Some(4), None    , 4..12),
        case::front_more_back_none(12, 4, Some(5), None    , 5..12),
        case::front_none_back_part(12, 4, None   , Some(11), 0..11),
        case::front_none_back_full(12, 4, None   , Some( 8), 0.. 8),
        case::front_none_back_more(12, 4, None   , Some( 7), 0.. 7),
        case::front_part_back_part(12, 4, Some(2), Some( 9), 2.. 9),
        case::front_part_back_full(12, 4, Some(2), Some( 8), 2.. 8),
        case::front_part_back_more(12, 4, Some(2), Some( 6), 2.. 6),
        case::front_full_back_part(12, 4, Some(4), Some( 9), 4.. 9),
        case::front_full_back_full(12, 4, Some(4), Some( 8), 4.. 8),
        case::front_full_back_more(12, 4, Some(4), Some( 6), 4.. 6),
        case::front_more_back_part(12, 4, Some(5), Some( 9), 5.. 9),
        case::front_more_back_full(12, 4, Some(5), Some( 8), 5.. 8),
        case::front_more_back_more(12, 4, Some(5), Some( 6), 5.. 6),
    )]
    fn test_name(
        #[case] d: usize,
        #[case] c: usize,
        #[case] min: Option<usize>,
        #[case] max: Option<usize>,
        #[case] expected: std::ops::Range<usize>,
    ) {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("chunking-test.h5");
        let dataset = "stuff";
        write_chunked(d, c, &file_path).unwrap();
        let read = iter_dataset::<f32>(&file_path, dataset, Bounds { min, max })
            .unwrap()
            .collect_vec();
        let expected = expected.map(|n| n as f32).collect_vec();
        println!("{read:?}");
        println!("{expected:?}");
        assert_eq!(read, expected);
    }

    // Only used for testing, but maybe should be expanded to a more general tool
    fn write_chunked(data_size: usize, chunk_size: usize, filename: impl AsRef<Path>) -> hdf5::Result<()> {
        let file = hdf5::File::create(filename)?;

        let ds = file
            .new_dataset::<f32>()
            .chunk(chunk_size)
            .shape(0..)
            .create("stuff")?;

        let data = (0..data_size).map(|i| i as f32);
        for chunk in &data.chunks(chunk_size) {
            let chunk = chunk.collect_vec();
            let old_size = ds.shape()[0];
            ds.resize(old_size + chunk_size)?;
            ds.write_slice(&chunk, old_size..)?
        }

        Ok(())
    }

}

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
        let read_data = read_dataset::<Outer>(&file_path, "just-testing/nested", Bounds::none())?.to_vec();
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
