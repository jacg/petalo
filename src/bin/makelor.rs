use structopt::StructOpt;
use itertools::Itertools;
use indicatif::{ProgressBar, ProgressStyle};
use petalo::{io, weights::LOR};
use petalo::io::hdf5::{SensorXYZ, Hdf5Lor};
use petalo::types::{Point, Time};

#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "makelor", about = "Create LORs from MC data")]
pub struct Cli {
    /// HDF5 input files with waveform and charge tables
    pub infiles: Vec<String>,

    /// HDF5 output file for LORs found in input file
    #[structopt(short, long)]
    pub out: String,

    /// Ignore sensors which do not receive this number of photons
    #[structopt(short, long, default_value = "4")]
    pub threshold: u64,

    /// Print LORs on stdout
    #[structopt(short, long)]
    pub print: bool,

    // TODO allow using different group/dataset in output
}

fn main() -> hdf5::Result<()> {
    let args = Cli::from_args();
    // --- Progress bar --------------------------------------------------------------
    let files_pb = ProgressBar::new(args.infiles.len() as u64).with_message(args.infiles[0].clone());
    files_pb.set_style(ProgressStyle::default_bar()
                       .template("Reading file: {msg}\n[{elapsed_precise}] {wide_bar} {pos}/{len} ({eta})")
        );
    files_pb.tick();
    // --- Process input files -------------------------------------------------------
    let mut xyzs = None; // sensor positions
    let mut lors: Vec<Hdf5Lor> = vec![];
    let threshold = args.threshold;
    let mut n_events = 0;
    let mut failed_files = vec![];
    for infile in args.infiles {
        files_pb.set_message(infile.clone());
        if let Ok(qts) = read_file(&infile, &mut xyzs) {
            let qts = qts.into_iter().filter(|h| h.q > threshold).collect(); // TODO: pass iterator to group_by_event
            let events = group_by_event(qts);
            for hits in events {
                n_events += 1;
                if let Some(lor) = lor(&hits, xyzs.as_ref().unwrap()) {
                    lors.push(lor.into());
                }
            }
        } else { failed_files.push(infile); }
        files_pb.inc(1);
    }
    files_pb.finish_with_message("<finished reading files>");
    println!("{} / {} ({}%) events produced LORs", lors.len(), n_events,
             100 * lors.len() / n_events);
    // --- write lors to hdf5 --------------------------------------------------------
    println!("Writing LORs to {}", args.out);
    hdf5::File::create(args.out)?
        .create_group("reco_info")?
        .new_dataset::<Hdf5Lor>().create("lors", lors.len())?
        .write(&lors)?;
    // --- Report any files that failed no be read -----------------------------------
    if !failed_files.is_empty() {
        println!("Warning: failed to read:");
        for file in failed_files {
            println!("  {}", file);
        }
    }
    Ok(())
}

fn lor(hits: &[QT0], xyzs: &SensorXyz) -> Option<LOR> {
    let (cluster_a, cluster_b) = group_into_clusters(hits, xyzs)?;
    //println!("{} + {} = {} ", cluster_a.len(), cluster_b.len(), hits.len());
    let (pa, ta) = cluster_xyzt(&cluster_a, &xyzs)?;
    let (pb, tb) = cluster_xyzt(&cluster_b, &xyzs)?;
    //println!("{:?} {:?}", xyzt_a, xyzt_b);
    Some(LOR::new(ta, tb, pa, pb))
}

fn cluster_xyzt(hits: &[QT0], xyzs: &SensorXyz) -> Option<(Point, Time)> {
    let (x,y,z) = barycentre(hits, xyzs)?;
    let ts = k_smallest(10, hits.into_iter().map(|h| h.t0))?;
    let t = mean(&ts)?;
    //let t = hits.iter().cloned().map(|h| Finite::<f32>::from(h.t0)).min()?.into();
    Some((Point::new(x,y,z),t))
}

fn mean(data: &Vec<f32>) -> Option<f32> {
    let n = data.len();
    if n > 0 {
        let sum: f32 = data.iter().sum();
        Some(sum / n as f32)
    }
    else { None }
}

// TODO: I'm sure it can be done more efficiently than this quick hack
fn k_smallest<I>(k: usize, data: I) -> Option<Vec<f32>>
where
    I: IntoIterator<Item = f32>,
{
    use ordered_float::NotNan;
    let mut heap = std::collections::BinaryHeap::<NotNan<f32>>::new();
    for v in data { heap.push(NotNan::new(-v).ok()?); } // NB negate to turn max-heap into min-heap
    let mut result = vec![];
    for _ in 0..k {
        let x = heap.pop()?;
        result.push((-x).into()); // NB undo earlier negation
    }
    Some(result)
}

fn find_sensor_with_highest_charge(sensors: &[QT0]) -> Option<u64> {
    sensors.iter().max_by_key(|e| e.q).map(|e| e.sensor_id)
}

type SensorXyz = std::collections::HashMap<u64, (f32, f32, f32)>;

fn group_into_clusters(hits: &[QT0], xyzs: &SensorXyz) -> Option<(Vec::<QT0>, Vec::<QT0>)> {
    let sensor_with_highest_charge = find_sensor_with_highest_charge(hits)?;
    let mut a = Vec::<QT0>::new();
    let mut b = Vec::<QT0>::new();
    let &(xm, ym, _) = xyzs.get(&sensor_with_highest_charge)?;
    for hit in hits.iter().cloned() {
        let &(x, y, _) = xyzs.get(&hit.sensor_id)?;
        if dot((xm, ym), (x,y)) > 0.0 { a.push(hit) }
        else                          { b.push(hit) }
    }
    if b.len() > 0 { Some((a, b)) } else { None }
}

fn dot((x1,y1): (f32, f32), (x2,y2): (f32, f32)) -> f32 { x1*x2 + y1*y2 }

fn barycentre(hits: &[QT0], xyzs: &SensorXyz) -> Option<(f32, f32, f32)> {
    if hits.len() == 0 { return None }
    let mut qs = 0_f32;
    let mut xx = 0.0;
    let mut yy = 0.0;
    let mut zz = 0.0;
    for QT0{ sensor_id, q, .. } in hits {
        let (x, y, z) = xyzs.get(&sensor_id)?;
        let q = *q as f32;
        qs += q;
        xx += x * q;
        yy += y * q;
        zz += z * q;
    }
    Some((xx / qs, yy / qs, zz / qs))
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct QT0 {
    pub event_id: u64,
    pub sensor_id: u64,
    pub q: u64,
    pub t0: f32,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Waveform {
    pub event_id: u64,
    pub sensor_id: u64,
    pub time: f32,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Qtot {
    pub event_id: u64,
    pub sensor_id: u64,
    pub charge: u64,
}

type SensorMap = std::collections::HashMap<u64, (f32, f32, f32)>;

// Nasty output parameter inteface ...
fn read_file(infile: &String, xyzs: &mut Option<SensorMap>) -> hdf5::Result<Vec<QT0>> {
    // Read sensor configuration, if not yet set.
    if xyzs.is_none() {
        // TODO: We're assuming that that the sensor configurations in all files
        // are the same. We should verify it.
        let yy = io::hdf5::read_table::<SensorXYZ>(&infile, "MC/sensor_xyz"  , None)?;
        let mut xx = vec![];
        xx.extend_from_slice(yy.as_slice().unwrap());
        // TODO ndarray 0.14 -> 0.15: breaks our code in hdf5
        // joined.extend_from_slice(data.into_slice());
        *xyzs = Some(make_sensor_position_map(xx));
    }
    // Read charges and waveforms
    let qs = io::hdf5::read_table::<Qtot     >(&infile, "MC/total_charge", None)?;
    let ts = io::hdf5::read_table::<Waveform >(&infile, "MC/waveform"    , None)?;
    Ok(combine_tables(qs, ts))
}

fn combine_tables(qs: ndarray::Array1<Qtot>, ts: ndarray::Array1<Waveform>) -> Vec<QT0> {
    let mut qts = vec![];
    let mut titer = ts.into_iter();
    for &Qtot{ event_id, sensor_id, charge:q} in qs.iter() {
        while let Some(&Waveform{ event_id: te, sensor_id: ts, time:t}) = titer.next() {
            if event_id == te && sensor_id == ts {
                qts.push(QT0{ event_id, sensor_id, q, t0:t });
                break;
            }
        }
    }
    qts
}

fn group_by_event(qts: Vec<QT0>) -> Vec<Vec<QT0>> {
    qts.into_iter()
        .group_by(|h| h.event_id)
        .into_iter()
        .map(|(_, group)| group.collect())
        .collect()
}

fn make_sensor_position_map(xyzs: Vec<SensorXYZ>) -> SensorMap {
    xyzs.iter().cloned()
        .map(|SensorXYZ{sensor_id, x, y, z}| (sensor_id, (x as f32, y as f32, z as f32)))
        .collect()
}
