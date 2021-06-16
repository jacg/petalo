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
    pub out: Option<String>,

    /// Ignore sensors which do not receive this number of photons
    #[structopt(short, long, default_value = "4")]
    pub threshold: u64,

    /// Print LORs on stdout
    #[structopt(short, long)]
    pub print: bool,

    // TODO allow using different group/dataset in output
}

// Nasty output parameter inteface ...
fn read_file(infile: &String, xyzs: &mut Option<Vec<SensorXYZ>>) -> hdf5::Result<Vec<QT0>> {
    //println!("Reading file {}", infile);
    // Read sensor configuration, if not yet set.
    if xyzs.is_none() {
        // TODO: We're assuming that that the sensor configurations in all files
        // are the same. We should verify it.
        let yy = io::hdf5::read_table::<SensorXYZ>(&infile, "MC/sensor_xyz"  , None)?;
        let mut xx = vec![];
        xx.extend_from_slice(yy.as_slice().unwrap());
        // TODO ndarray 0.14 -> 0.15: breaks our code in hdf5
        // joined.extend_from_slice(data.into_slice());
        *xyzs = Some(xx);
    }
    // Read charges and waveforms
    let qs = io::hdf5::read_table::<Qtot     >(&infile, "MC/total_charge", None)?;
    let ts = io::hdf5::read_table::<Waveform >(&infile, "MC/waveform"    , None)?;
    Ok(combine_tables(qs, ts))
}

fn combine_tables(qs: ndarray::Array1<Qtot>, ts: ndarray::Array1<Waveform>) -> Vec<QT0> {
    //println!("Combining tables ... ");
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
    //println!("done");
    qts
}

fn group_by_event(qts: Vec<QT0>) -> Vec<Vec<QT0>> {
    qts.into_iter()
        .group_by(|h| h.event_id)
        .into_iter()
        .map(|(_, group)| group.collect())
        .collect()
}

fn main() -> hdf5::Result<()> {
    let args = Cli::from_args();
    // Progress bars
    let files_pb = ProgressBar::new(args.infiles.len() as u64)
        .with_message(args.infiles[0].clone());
    files_pb.set_style(ProgressStyle::default_bar()
                       .template("Reading file: {msg}\n[{elapsed_precise}] {wide_bar} {pos}/{len} ({eta})")
        );
    files_pb.tick();
    // --- read data -----------------------------------------------------------------
    let mut xyzs = None;
    let mut qts  = vec![];
    let threshold = args.threshold;
    for infile in args.infiles {
        files_pb.set_message(infile.clone());
        qts.extend(read_file(&infile, &mut xyzs)?.into_iter()
                   // ignore sensors with small number of hits
                   .filter(|h| h.q > threshold)
        );
        files_pb.inc(1);
    }
    files_pb.finish_with_message("");
    // --- make map of sensor x-y positions ------------------------------------------
    let xyzs = xyzs.unwrap().iter().cloned()
        .map(|SensorXYZ{sensor_id, x, y, z}| (sensor_id, (x as f32, y as f32, z as f32)))
        .collect::<std::collections::HashMap<_,_>>();
    // --- group data into events ----------------------------------------------------
    let events: Vec<Vec<QT0>> = group_by_event(qts);
    let n_events = events.len();
    // --- calculate lors for each event ---------------------------------------------
    println!("Calculating lors");
    let lors_pb = ProgressBar::new(events.len() as u64);
    let mut count_interesting_events = 0_usize;
    let mut lors = Vec::<Hdf5Lor>::new();
    let write = args.out.is_some();
    for hits in events {
        if let Some(lor) = lor(&hits, &xyzs) {
            count_interesting_events += 1;
            if args.print { println!("{:6} {}", count_interesting_events-1, lor) };
            if      write { lors.push(lor.into()) };
        }
        lors_pb.inc(1);
    }
    lors_pb.finish_with_message("done");
    println!("{} / {} ({}%) have 2 clusters", count_interesting_events, n_events,
             100 * count_interesting_events / n_events);
    // --- write lors to hdf5 -------------------------------------------------------
    if let Some(file_name) = args.out {
        println!("Writing LORs to {}", file_name);
        hdf5::File::create(file_name)?
            .create_group("reco_info")?
            .new_dataset::<Hdf5Lor>().create("lors", lors.len())?
            .write(&lors)?;
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
