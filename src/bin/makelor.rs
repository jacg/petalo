use structopt::StructOpt;
use petalo::{io, weights::LOR};
use petalo::io::hdf5::{SensorXYZ, Hdf5Lor};
use petalo::types::{Point, Time};

#[derive(StructOpt, Debug, Clone)]
#[structopt(name = "makelor", about = "Create LORs from MC data")]
pub struct Cli {
    /// HDF5 input file with (event_id, sensor_id, Q, t0)
    pub infile: String,

    /// HDF5 output file for LORs found in input file
    #[structopt(short, long)]
    pub out: Option<String>,

    /// Ignore sensors which do not receive this number of photons
    #[structopt(short, long, default_value = "4")]
    pub threshold: u64,

    /// Print LORs on stdout
    #[structopt(short, long)]
    pub print: bool,

}

fn main() -> hdf5::Result<()> {
    let args = Cli::from_args();
    // --- read data -----------------------------------------------------------------
    let infile = args.infile;
    let xyzs = io::hdf5::read_table::<SensorXYZ>(&infile, "MC/sensor_xyz", None)?;
    let qts  = io::hdf5::read_table::<QT0      >(&infile, "MC/q_t0"      , None)?;
    // --- ignore sensors with small number of hits ----------------------------------
    let threshold = args.threshold;
    let qts = qts.iter().filter(|h| h.q > threshold).cloned().collect::<Vec<_>>();

    // --- make map of sensor x-y positions ------------------------------------------
    let xyzs = xyzs.iter().cloned()
        .map(|SensorXYZ{sensor_id, x, y, z}| (sensor_id, (x as f32, y as f32, z as f32)))
        .collect::<std::collections::HashMap<_,_>>();
    // --- find available event numbers ----------------------------------------------
    let event_numbers = qts.iter().map(|e| e.event_id).collect::<std::collections::BTreeSet<_>>();
    // --- calculate lors for each event ---------------------------------------------
    let mut count_interesting_events = 0_u16;
    let mut lors = Vec::<Hdf5Lor>::new();
    let write = args.out.is_some();
    for e in event_numbers.iter().cloned() { // TODO replace filter with groupby
        let hits = qts.iter().filter(|h| h.event_id == e).cloned().collect::<Vec<_>>();
        if let Some(lor) = lor(&hits, &xyzs) {
            count_interesting_events += 1;
            if args.print { println!("{:6} {}", count_interesting_events-1, lor) };
            if      write { lors.push(lor.into()) };
        }
    }
    println!("{} / {} have 2 clusters", count_interesting_events, event_numbers.len());
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