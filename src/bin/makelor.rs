use structopt::StructOpt;
use itertools::Itertools;
use indicatif::{ProgressBar, ProgressStyle};
use petalo::io;
use petalo::io::hdf5::{SensorXYZ, Hdf5Lor};
use petalo::Energyf32;
use petalo::{Length, Time, Point, Ratio};
use geometry::uom::mmps::f32::Area;
use geometry::uom::uomcrate::ConstZero;
use petalo::utils::group_digits;

// TODO: try to remove the need for these
use geometry::uom::{mm, mm_, ns, ns_, ratio};
// The mair problems seems to be that uom types do not implement various third-party traits, such as:
// + std::iter::Sum
// + hdf5:H5Type
// + From<ordered_float::NotNan>
// + Dbscan something or other

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "makelor", about = "Create LORs from MC data")]
pub struct Cli {
    /// HDF5 input files with waveform and charge tables
    pub infiles: Vec<String>,

    /// HDF5 output file for LORs found in input file
    #[structopt(short, long)]
    pub out: String,

    #[structopt(subcommand)]
    reco: Reco,

    // TODO allow using different group/dataset in output
}

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
enum Reco {

    /// Reconstruct LORs from first vertices in LXe
    #[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
    FirstVertex,

    /// Reconstruct LORs from barycentre of vertices in LXe
    #[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
    BaryVertex,

    /// Reconstruct LORs from clusters found by splitting cylinder in half
    #[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
    Half {
        /// Ignore sensors with fewer hits
        #[structopt(short, long = "charge-threshold", default_value = "4")]
        q: u32,
    },

    /// Reconstruct LORs form DBSCAN clusters
    #[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
    Dbscan {
        /// Ignore sensors with fewer hits
        #[structopt(short, long = "charge-threshold", default_value = "4")]
        q: u32,

        /// Minimum number of sensors in cluster
        #[structopt(short = "n", long, default_value = "10")]
        min_count: usize,

        /// Maximum distance between neighbours in cluster
        #[structopt(short = "d", long, default_value = "100 mm")]
        max_distance: Length,
    }
}

fn main() -> hdf5::Result<()> {
    let args = Cli::from_args();

    // Before starting the potentially long computation, make sure that we can
    // write the result to the requested destination. If the directory where
    // results will be written does not exist yet, make it. Panic if that's impossible.
    std::fs::create_dir_all(std::path::PathBuf::from(&args.out).parent().unwrap())
        .unwrap_or_else(|e| panic!("\n\nCan't write to \n\n   {}\n\n because \n\n   {e}\n\n", args.out));
    // --- Progress bar --------------------------------------------------------------
    let files_pb = ProgressBar::new(args.infiles.len() as u64).with_message(args.infiles[0].clone());
    files_pb.set_style(ProgressStyle::default_bar()
                       .template("Processing file: {msg}\n[{elapsed_precise}] {wide_bar} {pos}/{len} ({eta_precise})")
        );
    files_pb.tick();
    // --- Process input files -------------------------------------------------------
    let xyzs = read_sensor_map(&args.infiles[0])?;
    let mut lors: Vec<Hdf5Lor> = vec![];
    let mut n_events = 0;
    let mut failed_files = vec![];

    type FilenameToLorsFunction = Box<dyn Fn(&String) -> hdf5::Result<(Vec<Hdf5Lor>, usize)>>;
    let makelors: FilenameToLorsFunction = match args.reco {
        Reco::FirstVertex => Box::new(
            |infile: &String| -> hdf5::Result<(Vec<Hdf5Lor>, usize)> {
                let vertices = read_vertices(infile)?;
                let events = group_by(|v| v.event_id, vertices.into_iter());
                Ok((lors_from(&events, lor_from_first_vertices), events.len()))
            }),

        Reco::BaryVertex => Box::new(
            |infile: &String| -> hdf5::Result<(Vec<Hdf5Lor>, usize)> {
                let vertices = read_vertices(infile)?;
                let events = group_by(|v| v.event_id, vertices.into_iter());
                Ok((lors_from(&events, lor_from_barycentre_of_vertices), events.len()))
            }),

        Reco::Half{q} => Box::new(
            move |infile: &String| -> hdf5::Result<(Vec<Hdf5Lor>, usize)> {
                let qts = read_qts(infile)?;
                let events = group_by(|h| h.event_id, qts.into_iter().filter(|h| h.q >= q));
                Ok((lors_from(&events, |evs| lor_from_hits(evs, &xyzs)), events.len()))
            }),

        Reco::Dbscan { q, min_count, max_distance } => Box::new(
            move |infile: &String| -> hdf5::Result<(Vec<Hdf5Lor>, usize)> {
                let qts = read_qts(infile)?;
                let events = group_by(|h| h.event_id, qts.into_iter().filter(|h| h.q >= q));
                Ok((lors_from(&events, |evs| lor_from_hits_dbscan(evs, &xyzs, min_count, max_distance)), events.len()))
            }),
    };


    for infile in args.infiles {
        // TODO message doesn't appear until end of iteration
        files_pb.set_message(format!("{}. Found {} LORs in {} events, so far ({}%).",
                                     infile.clone(), group_digits(lors.len()), group_digits(n_events),
                                     if n_events > 0 {100 * lors.len() / n_events} else {0}));
        if let Ok((new_lors, envlen)) = makelors(&infile) {
            n_events += envlen;
            lors.extend_from_slice(&new_lors);
        } else { failed_files.push(infile); }
        files_pb.inc(1);

    }
    files_pb.finish_with_message("<finished processing files>");
    println!("{} / {} ({}%) events produced LORs", group_digits(lors.len()), group_digits(n_events),
             100 * lors.len() / n_events);
    // --- write lors to hdf5 --------------------------------------------------------
    println!("Writing LORs to {}", args.out);
    hdf5::File::create(args.out)?
        .create_group("reco_info")?
        .new_dataset_builder()
        .with_data(&lors)
        .create("lors")?;
    // --- Report any files that failed no be read -----------------------------------
    if !failed_files.is_empty() {
        println!("Warning: failed to read the following files:");
        for file in failed_files.iter() {
            println!("  {}", file);
        }
        let n = failed_files.len();
        let plural = if n == 1 { "" } else { "s" };
        println!("Warning: failed to read {} file{}:", n, plural);
    }
    Ok(())
}

fn lors_from<T>(events: &[Vec<T>], mut one_lor: impl FnMut(&[T]) -> Option<Hdf5Lor>) -> Vec<Hdf5Lor> {
    events.iter()
        .flat_map(|data| one_lor(data))
        .collect()
}

#[allow(nonstandard_style)]
fn lor_from_first_vertices(vertices: &[Vertex]) -> Option<Hdf5Lor> {
    let mut in_lxe = vertices.iter().filter(|v| v.volume_id == 0);
    let &Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: E2, ..} = in_lxe.find(|v| v.track_id == 2)?;
    let &Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: E1, ..} = in_lxe.find(|v| v.track_id == 1)?;
    Some(Hdf5Lor {
        dt: t2 - t1,                   x1, y1, z1,   x2, y2, z2,
        q1: f32::NAN, q2: f32::NAN,        E1,           E2,
    })
}

#[allow(nonstandard_style)]
fn lor_from_barycentre_of_vertices(vertices: &[Vertex]) -> Option<Hdf5Lor> {
    let (a,b): (Vec<_>, Vec<_>) = vertices
        .iter()
        .filter(|v| v.volume_id == 0)
        .partition(|v| v.track_id == 1);

    let Barycentre { x: x1, y: y1, z: z1, t: t1, E: E1 } = vertex_barycentre(&a)?;
    let Barycentre { x: x2, y: y2, z: z2, t: t2, E: E2 } = vertex_barycentre(&b)?;

    let nan = f32::NAN; let q1 = nan; let q2 = nan;

    Some(Hdf5Lor {
        dt: ns_(t2 - t1),
        x1: mm_(x1), y1: mm_(y1), z1: mm_(z1),
        x2: mm_(x2), y2: mm_(y2), z2: mm_(z2),
        q1, q2, E1, E2,
    })
}

fn lor_from_hits(hits: &[QT], xyzs: &SensorMap) -> Option<Hdf5Lor> {
    let (cluster_a, cluster_b) = group_into_clusters(hits, xyzs)?;
    //println!("{} + {} = {} ", cluster_a.len(), cluster_b.len(), hits.len());
    let (p1, t1) = cluster_xyzt(&cluster_a, xyzs)?;
    let (p2, t2) = cluster_xyzt(&cluster_b, xyzs)?;
    //println!("{:?} {:?}", xyzt_a, xyzt_b);
    Some(Hdf5Lor {
        dt: ns_(t2 - t1),
        x1: mm_(p1.x), y1: mm_(p1.y), z1: mm_(p1.z),
        x2: mm_(p2.x), y2: mm_(p2.y), z2: mm_(p2.z),
        // TODO qs and Es missing
        q1: f32::NAN, q2: f32::NAN,   E1: f32::NAN, E2: f32::NAN,
    })
}

fn xxx(labels: &ndarray::Array1<Option<usize>>) -> usize {
    labels.iter()
        .flat_map(|l| *l)
        .max()
        .map_or(0, |n| n+1)
}

#[cfg(test)]
mod test_n_clusters {
    use super::*;
    use ndarray::array;
    #[test]
    fn test_n_clusters() {
        let a = array![None, None];
        assert_eq!(xxx(&a), 0);
        let _a = array![Some(2)];
    }
}

fn lor_from_hits_dbscan(hits: &[QT], xyzs: &SensorMap, min_points: usize, tolerance: Length) -> Option<Hdf5Lor> {
    use linfa_clustering::AppxDbscan;
    use linfa::traits::Transformer;
    let active_sensor_positions: ndarray::Array2<f32> = hits.iter()
        .flat_map(|QT { sensor_id, ..}| xyzs.get(sensor_id))
        .map(|&(x,y,z)| [mm_(x), mm_(y), mm_(z)])
        .collect::<Vec<_>>()
        .into();
    let params = AppxDbscan::params(min_points)
        .tolerance(mm_(tolerance)); // > 7mm between sipm centres
    let labels = params.transform(&active_sensor_positions).ok()?;
    let n_clusters = xxx(&labels);
    if n_clusters != 2 { return None }
    let mut cluster: [Vec<f32>; 2] = [vec![], vec![]];
    for (c, point) in labels.iter().zip(active_sensor_positions.outer_iter()) {
        if let Some(c) = c {
            cluster[*c].extend(&point);
        }
    }
    fn cluster_centroid(vec: Vec<f32>) -> Option<ndarray::Array1<f32>> {
        let it = ndarray::Array2::from_shape_vec((vec.len()/3, 3), vec.clone()).unwrap()
            .mean_axis(ndarray::Axis(0));
        if it.is_none() {
            println!("Failed to find centroid of cluster {:?}", vec);
        }
        it
    }
    let a = cluster_centroid(cluster[0].clone())?;
    let b = cluster_centroid(cluster[1].clone())?;
    Some(Hdf5Lor {
        dt: 0.0, // TODO dt missing
        x1: a[0], y1: a[1], z1: a[2],
        x2: b[0], y2: b[1], z2: b[2],
        // TODO qs and Es missing
        q1: f32::NAN, q2: f32::NAN,   E1: f32::NAN, E2: f32::NAN,
    })
}

fn cluster_xyzt(hits: &[QT], xyzs: &SensorMap) -> Option<(Point, Time)> {
    let (x,y,z) = sipm_charge_barycentre(hits, xyzs)?;
    let ts = k_smallest(10, hits.iter().map(|h| h.t))?;
    let t = mean(&ts)?;
    //let t = hits.iter().cloned().map(|h| Finite::<f32>::from(h.t0)).min()?.into();
    Some((Point::new(x,y,z),t))
}

fn mean(data: &[Time]) -> Option<Time> {
    let n = data.len();
    if n > 0 {
        let sum: Time = data.iter().cloned().sum();
        Some(sum / n as f32)
    }
    else { None }
}

// TODO: I'm sure it can be done more efficiently than this quick hack
fn k_smallest<I>(k: usize, data: I) -> Option<Vec<Time>>
where
    I: IntoIterator<Item = Time>,
{
    use ordered_float::NotNan;
    let mut heap = std::collections::BinaryHeap::<NotNan<f32>>::new();
    for v in data { heap.push(NotNan::new(-ns_(v)).ok()?); } // NB negate to turn max-heap into min-heap
    let mut result = vec![];
    for _ in 0..k {
        let x = heap.pop()?;
        result.push(ns((-x).into())); // NB undo earlier negation
    }
    Some(result)
}

fn find_sensor_with_highest_charge(sensors: &[QT]) -> Option<u32> {
    sensors.iter().max_by_key(|e| e.q).map(|e| e.sensor_id)
}

type SensorMap = std::collections::HashMap<u32, (Length, Length, Length)>;

fn group_into_clusters(hits: &[QT], xyzs: &SensorMap) -> Option<(Vec::<QT>, Vec::<QT>)> {
    let sensor_with_highest_charge = find_sensor_with_highest_charge(hits)?;
    let mut a = Vec::<QT>::new();
    let mut b = Vec::<QT>::new();
    let &(xm, ym, _) = xyzs.get(&sensor_with_highest_charge)?;
    for hit in hits.iter().cloned() {
        let &(x, y, _) = xyzs.get(&hit.sensor_id)?;
        if dot((xm, ym), (x,y)) > Area::ZERO { a.push(hit) }
        else                                 { b.push(hit) }
    }
    if !b.is_empty() { Some((a, b)) } else { None }
}

fn dot((x1,y1): (Length, Length), (x2,y2): (Length, Length)) -> Area { x1*x2 + y1*y2 }

fn sipm_charge_barycentre(hits: &[QT], xyzs: &SensorMap) -> Option<(Length, Length, Length)> {
    if hits.is_empty() { return None }
    let mut qs = Ratio::ZERO;
    let mut xx = Length::ZERO;
    let mut yy = Length::ZERO;
    let mut zz = Length::ZERO;
    for &QT{ sensor_id, q, .. } in hits {
        let &(x, y, z) = xyzs.get(&sensor_id)?;
        let q = ratio(q as f32);
        qs += q;
        xx += x * q;
        yy += y * q;
        zz += z * q;
    }
    Some((xx / qs, yy / qs, zz / qs))
}

#[allow(nonstandard_style)]
#[derive(Copy, Clone)]
struct Barycentre {
    x: Length,
    y: Length,
    z: Length,
    t: Time,
    E: Energyf32,
}

#[allow(nonstandard_style)]
fn vertex_barycentre(vertices: &[&Vertex]) -> Option<Barycentre> {
    if vertices.is_empty() { return None }
    let mut delta_E  = 0_f32;
    let mut rr = Length::ZERO;
    let mut cc = Ratio::ZERO;
    let mut zz = Length::ZERO;
    let mut tt = Time::ZERO;
    for &&Vertex { x, y, z, t, pre_KE, post_KE, .. } in vertices {
        let dE = pre_KE - post_KE;
        let (x, y, z, t) = (mm(x), mm(y), mm(z), ns(t));
        delta_E += dE;
        let radius = (x*x + y*y).sqrt();
        rr += radius * dE;
        cc += (x / radius) * dE;
        zz += z * dE;
        tt += t * dE; // TODO figure out what *really* needs to be done here
    };
    rr /= delta_E;
    cc /= delta_E;
    let xx = rr * cc;
    let yy = rr * (ratio(1.0) - cc*cc).sqrt();
    Some(Barycentre { x: xx, y: yy, z: zz / delta_E, t: tt / delta_E, E: delta_E })
}

#[cfg(test)]
mod test_vertex_barycentre {
    use super::*;
    use float_eq::assert_float_eq;
    use geometry::uom::radian;
    use std::f32::consts::PI;

    /// Create a vertex with interesting x,y and dummy values elsewhere
    fn vertex(x: Length, y: Length) -> Vertex {
        Vertex {
            // interesting values
            x: mm_(x), y: mm_(y),
            // dummy values
            z: 23.4, t: 123.0,
            event_id: 0, parent_id: 0, track_id: 0, process_id: 0, volume_id: 0,
            moved: 0.0, deposited: 0, pre_KE: 511.0, post_KE: 0.0,
        }
    }

    #[test]
    fn curvature_should_not_reduce_radius() {

        // Place all vertices at same distance from axis
        let r = mm(355.0);

        // Distribute vertices on quarter circle
        let vertices: Vec<Vertex> = (0..10)
            .map(|i| radian(i as f32 * PI / 20.0))
            .map(|angle| (r * angle.cos(), r * angle.sin()))
            .map(|(x,y)| vertex(x,y))
            .collect();

        // Create vector of vertex refs, as required by vertex_barycentre
        let vertex_refs: Vec<&Vertex> = vertices.iter().collect();

        let Barycentre { x, y, .. } = vertex_barycentre(&vertex_refs).unwrap();
        let bary_r = (x*x + y*y).sqrt();
        assert_float_eq!(mm_(bary_r), mm_(r), ulps <= 1);
    }

    #[test]
    fn angular_distribution_mean() {

        let r = 100.0;

        // One vertex at 0 degrees, the other at 90 degrees.
        let vertices = vec![vertex(mm(100.0), mm(0.0)), vertex(mm(0.0), mm(100.0))];

        // The barycentre of the above points should be at 45 degrees or pi / 4.
        let angle = PI / 4.0;
        let (expected_x, expected_y) = (r * angle.cos(), r * angle.sin());

        // Create vector of vertex refs, as required by vertex_barycentre
        let vertex_refs: Vec<&Vertex> = vertices.iter().collect();

        let Barycentre { x, y, .. } = vertex_barycentre(&vertex_refs).unwrap();

        assert_float_eq!(mm_(x), expected_x, ulps <= 1);
        assert_float_eq!(mm_(y), expected_y, ulps <= 1);
    }
}

#[derive(Clone, Debug)]
pub struct QT {
    pub event_id: u32,
    pub sensor_id: u32,
    pub q: u32,
    pub t: Time,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
#[allow(nonstandard_style)]
pub struct Vertex {
    event_id: u32,
    track_id: u32,
    parent_id: u32,
    x: f32,
    y: f32,
    z: f32,
    t: f32,
    moved: f32,
    pre_KE: f32,
    post_KE: f32,
    deposited: u32,
    process_id: u32, // NB these may differ across
    volume_id: u32,  // different files
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Waveform {
    pub event_id: u32,
    pub sensor_id: u32,
    pub time: f32,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct Qtot {
    pub event_id: u32,
    pub sensor_id: u32,
    pub charge: u32,
}

// TODO Is there really no simpler way?
fn array_to_vec<T: Clone>(array: ndarray::Array1<T>) -> Vec<T> {
    let mut vec = vec![];
    vec.extend_from_slice(array.as_slice().unwrap());
    // TODO ndarray 0.14 -> 0.15: breaks our code in hdf5
    // joined.extend_from_slice(data.into_slice());
    vec
}



fn read_sensor_map(filename: &str) -> hdf5::Result<SensorMap> {
    // TODO: refactor and hide in a function
    let array = io::hdf5::read_table::<SensorXYZ>(filename, "MC/sensor_xyz"  , None)?;
    Ok(make_sensor_position_map(array_to_vec(array)))
}

fn read_vertices(filename: &str) -> hdf5::Result<Vec<Vertex>> {
    Ok(array_to_vec(io::hdf5::read_table::<Vertex>(filename, "MC/vertices", None)?))
}

fn read_qts(infile: &str) -> hdf5::Result<Vec<QT>> {
    // Read charges and waveforms
    let qs = io::hdf5::read_table::<Qtot     >(infile, "MC/total_charge", None)?;
    let ts = io::hdf5::read_table::<Waveform >(infile, "MC/waveform"    , None)?;
    Ok(combine_tables(qs, ts))
}

fn combine_tables(qs: ndarray::Array1<Qtot>, ts: ndarray::Array1<Waveform>) -> Vec<QT> {
    let mut qts = vec![];
    let mut titer = ts.iter();
    for &Qtot{ event_id, sensor_id, charge:q} in qs.iter() {
        for &Waveform{ event_id: te, sensor_id: ts, time:t} in titer.by_ref() {
            if event_id == te && sensor_id == ts {
                qts.push(QT{ event_id, sensor_id, q, t: ns(t) });
                break;
            }
        }
    }
    qts
}

fn group_by<T>(group_by: impl FnMut(&T) -> u32, qts: impl IntoIterator<Item = T>) -> Vec<Vec<T>> {
    qts.into_iter()
        .group_by(group_by)
        .into_iter()
        .map(|(_, group)| group.collect())
        .collect()
}

fn make_sensor_position_map(xyzs: Vec<SensorXYZ>) -> SensorMap {
    xyzs.iter().cloned()
        .map(|SensorXYZ{sensor_id, x, y, z}| (sensor_id, (mm(x), mm(y), mm(z))))
        .collect()
}

// Proof of concept: nested compound hdf5 types
#[cfg(test)]
mod test_nested_compound_hdf5 {

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
        let mut test_data = vec![];
        test_data.push(Outer{ id:0, r#true:Inner{ a: 123, b: 4.56 }, reco:Inner{ a: 789, b: 0.12} });
        test_data.push(Outer{ id:1, r#true:Inner{ a: 132, b: 45.6 }, reco:Inner{ a: 798, b: 10.2} });

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
        let read_data = petalo::io::hdf5::read_table::<Outer>(&file_path, "just-testing/nested", None)?;
        let read_data = super::array_to_vec(read_data);
        assert_eq!(test_data, read_data);
        println!("Test table written to {}", file_path);
        Ok(())
    }
}

// Proof of concept compound HDF5 containing array
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
        let mut test_data = vec![];
        test_data.push(Waveform::new(0, 8, vec![1.2, 3.4]));
        test_data.push(Waveform::new(9, 3, vec![5.6, 7.8, 9.0]));
        let filename = "just-testing.h5";
        hdf5::File::create(filename)?
            .create_group("foo")?
            .new_dataset_builder()
            .with_data(&test_data)
            .create("bar")?;
        Ok(())
    }
}

#[cfg(test)]
mod test_dbscan {
    use linfa_clustering::{Dbscan, AppxDbscan, generate_blobs};
    use linfa::traits::Transformer;
    use ndarray::{Axis, array};
    use ndarray_rand::rand::SeedableRng;
    use rand_isaac::Isaac64Rng;
    //use approx::assert_abs_diff_eq;

    #[test]
    fn test_appx () {
        // Our random number generator, seeded for reproducibility
        let seed = 42;
        let mut rng = Isaac64Rng::seed_from_u64(seed);

        // `expected_centroids` has shape `(n_centroids, n_features)`
        // i.e. three points in the 2-dimensional plane
        let expected_centroids = array![[0., 1.], [-10., 20.], [-1., 10.]];
        // Let's generate a synthetic dataset: three blobs of observations
        // (100 points each) centered around our `expected_centroids`
        let _observations = generate_blobs(100, &expected_centroids, &mut rng);

        let observations = array![[1.0, 1.0, 3.0],
                                  [1.1, 1.1, 3.0],
                                  [0.9, 0.9, 3.0],
                                  [0.9, 1.1, 3.0],
                                  [1.1, 0.9, 3.0],
                                  [1.0, 0.9, 3.0],
                                  [0.9, 1.0, 3.0],
                                  [1.0, 1.1, 3.0],
                                  [1.1, 1.0, 3.0],

                                  [2.0, 2.0, 3.0],
                                  [2.1, 2.1, 3.0],
                                  [1.9, 1.9, 3.0],
                                  [1.9, 2.1, 3.0],
                                  [2.1, 1.9, 3.0],
                                  [2.0, 1.9, 3.0],
                                  [1.9, 2.0, 3.0],
                                  [2.0, 2.1, 3.0],
                                  [2.1, 2.0, 3.0],

        ];

        // Let's configure and run our AppxDbscan algorithm
        // We use the builder pattern to specify the hyperparameters
        // `min_points` is the only mandatory parameter.
        // If you don't specify the others (e.g. `tolerance`, `slack`)
        // default values will be used.
        let min_points = 3;
        let params = AppxDbscan::params(min_points)
            .tolerance(1.0)
            .slack(1e-3);
        // Let's run the algorithm!
        let labels = params.transform(&observations);
        // Points are `None` if noise `Some(id)` if belonging to a cluster.
        println!("\n\nobs\n{:?}\n\n", observations);
        println!("\n\npar\n{:?}\n\n", params);
        println!("\n\nlab\n{:?}\n\n", labels);
    }

    #[test]
    fn test_full () {
        // Our random number generator, seeded for reproducibility
        let seed = 42;
        let mut rng = Isaac64Rng::seed_from_u64(seed);

        // `expected_centroids` has shape `(n_centroids, n_features)`
        // i.e. three points in the 2-dimensional plane
        let expected_centroids = array![[0., 1.], [-10., 20.], [-1., 10.]];
        // Let's generate a synthetic dataset: three blobs of observations
        // (100 points each) centered around our `expected_centroids`
        let observations: ndarray::Array2<f64> = generate_blobs(100, &expected_centroids, &mut rng);

        // Let's configure and run our DBSCAN algorithm
        // We use the builder pattern to specify the hyperparameters
        // `min_points` is the only mandatory parameter.
        // If you don't specify the others (e.g. `tolerance`)
        // default values will be used.
        let min_points = 3;
        let labels = Dbscan::params(min_points)
            .tolerance(1.0)
            .transform(&observations)
            .unwrap();
        // Points are `None` if noise `Some(id)` if belonging to a cluster.

        println!("\n\nobs\n{:?}", observations);
        println!("\n\nclu\n{:?}", labels);
        //println!("{:?}", params);

        let mut veccy = [vec![], vec![], vec![]];

        for (c, point) in labels.iter().zip(observations.outer_iter()) {
            if let Some(c) = c {
                let point: Vec<f64> = point.iter().copied().collect();
                veccy[*c].push(point)
            }
        }

        println!("\n{:4.1?}\n\n", veccy[0]);
        println!("\n{:4.1?}\n\n", veccy[1]);
        println!("\n{:4.1?}\n\n", veccy[2]);



        let mut arry:[Vec<f64>; 3] = [vec![], vec![], vec![]];

        for (c, point) in labels.iter().zip(observations.outer_iter()) {
            if let Some(c) = c {
                arry[*c].extend(&point);
            }
        }

        fn foo(vec: Vec<f64>, width: usize) -> ndarray::Array2<f64> {
            ndarray::Array2::from_shape_vec((vec.len()/width, width), vec).unwrap()
        }

        println!("aaaaa {:?}", foo(arry[0].clone(), 2).mean_axis(Axis(0)));


        // use ndarray::prelude::*;
        // println!("{:?}", observations.mean_axis(Axis(0)));

        // use ndarray::stack;
        // println!("\n{:4.1?}\n\n", stack![Axis(0), arry[0]].mean_axis((Axis(0))));
        // println!("\n{:4.1?}\n\n", stack![Axis(0), arry[1]]);
        // println!("\n{:4.1?}\n\n", arry[2]);

    }
}
