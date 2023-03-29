mod cli;
mod progress;

fn main() -> hdf5::Result<()> {
    let args = Cli::parse();

    // Before starting the potentially long computation, make sure that we can
    // write the result to the requested destination. If the directory where
    // results will be written does not exist yet, make it. Panic if that's impossible.
    std::fs::create_dir_all(PathBuf::from(&args.out).parent().unwrap())
        .unwrap_or_else(|e| panic!("\n\nCan't write to \n\n   {}\n\n because \n\n   {e}\n\n", args.out.display()));
    // --- Progress bar --------------------------------------------------------------
    let progress = progress::Progress::new(&args.infiles);
    // --- Process input files -------------------------------------------------------
    let pool = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build().unwrap();
    let xyzs = io::hdf5::sensors::read_sensor_map(&args.infiles[0])?;
    let files = args.infiles.clone();
    macro_rules! go { ($a:expr, $b:expr, $c:expr) => { compose_steps(files, &progress, $a, $b, $c) }; }
    let lors: Vec<_> = pool.install( || match &args.reco {
        d@Reco::Discrete { .. } => go!(read_vertices, group_vertices, lor_from_discretized_vertices(d)),
        Reco::FirstVertex       => go!(read_vertices, group_vertices, lor_from_first_vertices),
        Reco::BaryVertex        => go!(read_vertices, group_vertices, lor_from_barycentre_of_vertices),
        &Reco::Half { q }       => go!(read_qts     , group_qts(q)  , lor_from_hits(&xyzs)),
        &Reco::Dbscan { q, min_count, max_distance } =>
            go!(read_qts, group_qts(q), lor_from_hits_dbscan(&xyzs, min_count, max_distance)),
        Reco::SimpleRec { .. } => todo!(), // Complex process related to obsolete detector design
    });

    // --- Report any files that failed no be read -----------------------------------
    progress.final_report();
    // --- write lors to hdf5 --------------------------------------------------------
    println!("Writing LORs to {}", args.out.display());
    hdf5::File::create(&args.out)?
        .create_group("reco_info")? // TODO rethink all the HDF% group names
        .new_dataset_builder()
        .with_data(&lors)
        .create("lors")?;
    Ok(())
}

fn compose_steps<T, E, G, M>(
    files: Vec<PathBuf>,
    stats: &Progress,
    extract_rows_from_table: E,
    group_by_event         : G,
    make_one_lor           : M,
) -> Vec<Hdf5Lor>
where
    T: Send,
    E: Fn(&Path)  -> hdf5::Result<Vec<T>> + Send + Sync,
    G: Fn(Vec<T>) -> Vec<Vec<T>>          + Send + Sync,
    M: Fn(  &[T]) -> Option<Hdf5Lor>      + Send + Sync,
{
    files
        .into_iter()
        .par_bridge() // With only 2 threads, locking on every LOR is cheap: reconsider if we ever overcome HDF5 global lock
        .map(|file| extract_rows_from_table(&file))      .inspect(|result| stats.read_file_done(result))
        .map(|rows| rows.unwrap_or_else(|_| vec![]))
        .map(group_by_event)                             .inspect(|x| stats.grouped(x))
        .flat_map_iter(|x| x.into_iter()
                       .filter_map(|x| make_one_lor(&x))).inspect(|x| stats.lor(x))
        .collect()
}

fn group_vertices(vertices: Vec<Vertex>) -> Vec<Vec<Vertex>> {
    group_by(|v| v.event_id, vertices)
}

fn group_qts(q: u32) -> impl Fn(Vec<QT>) -> Vec<Vec<QT>> {
    move |qts| {
        qts.into_iter()
           .filter(|h| h.q >= q)
           .group_by(|h| h.event_id)
           .into_iter()
           .map(|(_, group)| group.collect())
           .collect()
    }
}

fn read_vertices(infile: &Path) -> hdf5::Result<Vec<Vertex>> { io::hdf5::mc::read_vertices(infile, Bounds::none()) }
fn read_qts     (infile: &Path) -> hdf5::Result<Vec<  QT  >> { io::hdf5::sensors::read_qts(infile, Bounds::none()) }

// NOTE only using first vertex, for now
#[allow(nonstandard_style)]
fn lor_from_discretized_vertices(d: &Reco) -> impl Fn(&[Vertex]) -> Option<Hdf5Lor> + Copy {
    let &Reco::Discrete { r_min, dr, dz, da } = d else {
        panic!("lor_from_discretized_vertices called with variant other than Reco::Discrete")
    };
    use petalo::discrete::Discretize;
    let adjust = Discretize { r_min, dr, dz, da }.centre_of_nearest_box_fn_f32();
    move |vertices| {
        let mut in_scint = vertices_in_scintillator(vertices);

        // Optimistic version: grabs all remaining KE in first vertex
        let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: E2, ..} = in_scint.find(|v| v.track_id == 2)?;
        let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: E1, ..} = in_scint.find(|v| v.track_id == 1)?;

        // // Pessimistic version: picks vertex with highest dKE
        // use ordered_float::NotNan;
        // let in_scint = &mut vertices_in_scintillator(vertices); // Working around filter consuming the iterator
        // let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: b1, post_KE: a1, ..} = in_scint.filter(|v| v.track_id == 2).max_by_key(|v| NotNan::new(v.pre_KE - v.post_KE).unwrap())?;
        // let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: b2, post_KE: a2, ..} = in_scint.filter(|v| v.track_id == 2).max_by_key(|v| NotNan::new(v.pre_KE - v.post_KE).unwrap())?;
        // let (E1, E2) = (b1 - a1, b2 - a2);

        // let (ox1, oy1, oz1, ox2, oy2, oz2) = (x1, y1, z1, x2, y2, z2);
        let (x1, y1, z1) = adjust(x1, y1, z1);
        let (x2, y2, z2) = adjust(x2, y2, z2);

        // println!();
        // let dp1 = dist((ox1, oy1, oz1), (x1, y1, z1));
        // let dp2 = dist((ox2, oy2, oz2), (x2, y2, z2));
        // let mean_z = (z2+z1) / 2.0;
        // if E1.min(E2) < 500.0 { return None }
        // println!("{ox1:6.1} {oy1:6.1} {oz1:6.1}   {ox2:6.1} {oy2:6.1} {oz2:6.1}     {dp1:4.1}  {dp2:4.1}   {mean_z:6.1}");
        // println!("{x1:6.1} {y1:6.1} {z1:6.1}   {x2:6.1} {y2:6.1} {z2:6.1}    {E1:5.1} {E2:5.1}");

        Some(Hdf5Lor { dt: t2 - t1, x1, y1, z1, x2, y2, z2, q1: f32::NAN, q2: f32::NAN, E1, E2 })
    }
}

// fn dist((x1,y1,z1): (f32,f32,f32), (x2,y2,z2): (f32,f32,f32)) -> f32 {
//     (x2-x1).hypot(y2-y1).hypot(z2-z1)
// }
// ----- Imports -----------------------------------------------------------------------------------------
use std::path::{Path, PathBuf};
use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;
use units::{
    Length, Time, Ratio,
    todo::Energyf32, Area,
    uom::ConstZero
};
use petalo::{
    Point,
    config::mlem::Bounds,
    io::{self,
         hdf5::{
             Hdf5Lor,
             sensors::{QT, SensorMap},
             mc::Vertex
         }
    },
};
use cli::{Cli, Reco};
use progress::Progress;

// TODO: try to remove the need for these
use units::{mm, mm_, ns, ns_, ratio};
// The mair problems seems to be that uom types do not implement various third-party traits, such as:
// + std::iter::Sum
// + hdf5:H5Type
// + From<ordered_float::NotNan>
// + Dbscan something or other

// ----- TESTS ------------------------------------------------------------------------------------------
#[cfg(test)]
mod test_discretize {
    use rstest::rstest;
    use float_eq::assert_float_eq;

    use std::f32::consts::SQRT_2 as ROOT2;
    use std::f32::consts::TAU;

    #[rstest]
    //              rmin  dr     dz    da       x-in    y-in   z-in     x-out   y-out  z-out
    // On x/y-axis
    #[case::x_pos(( 90.0, 20.0,  3.0,  TAU), ( 109.9,    0.2,  28.6), ( 100.0 ,   0.0,  30.0))]
    #[case::x_neg(( 90.0, 20.0,  4.0,  TAU), (- 91.3,    0.2,  33.9), (-100.0 ,   0.0,  32.0))]
    #[case::y_pos(( 90.0, 20.0,  5.0,  TAU), (   0.2,   95.3,  50.1), (   0.0 , 100.0,  50.0))]
    #[case::y_neg(( 90.0, 20.0,  5.5,  TAU), (   0.2, -109.9,  52.3), (   0.0 ,-100.0,  55.0))]
    // Check that quadrants are preserved
    #[case::quad1((  1.9,  0.2,  1.0, 0.01), ( ROOT2,  ROOT2, -23.4), ( ROOT2, ROOT2, -23.0))]
    #[case::quad2((  1.9,  0.2,  1.0, 0.01), (-ROOT2,  ROOT2, - 9.9), (-ROOT2, ROOT2, -10.0))]
    #[case::quad3((  1.9,  0.2,  1.0, 0.01), (-ROOT2, -ROOT2, - 9.9), (-ROOT2,-ROOT2, -10.0))]
    #[case::quad4((  1.9,  0.2,  1.0, 0.01), ( ROOT2, -ROOT2,  23.4), ( ROOT2,-ROOT2,  23.0))]
    // Near x-axis
    #[case::xish1((280.0, 30.0,  4.0, 2.0 ), ( 309.9,    0.0, -31.6), ( 295.0,   0.0, -32.0))]
    #[case::xish2((280.0, 30.0,  5.0, 2.0 ), ( 280.1,    2.8,  31.6), ( 295.0,   2.0,  30.0))]
    #[case::xish3((280.0, 30.0,  6.0, 1.5 ), ( 294.6,   -0.8, -31.6), ( 295.0,  -1.5, -30.0))]
    // Other cases are very fiddly to verify by hand, but we should add some, in principle
    fn test_nearest_centre_of_box(
        #[case] detector: (f32, f32, f32, f32),
        #[case] vertex  : (f32, f32, f32),
        #[case] expected: (f32, f32, f32),
    ) {
        let (rmin, dr, dz, da) = detector;
        let (x, y, z) = vertex;
        let (x, y, z) = petalo::discrete::Discretize
            ::from_f32s_in_mm(rmin, dr, dz, da).
            centre_of_nearest_box_fn_f32()
            (x,y,z);
        assert_float_eq!((x,y,z), expected, abs <= (0.01, 0.01, 0.01));
    }
}

fn vertices_in_scintillator(vertices: &[Vertex]) -> impl Iterator<Item = Vertex> + '_ {
    vertices.iter().cloned().filter(|v| v.volume_id == 0)
}

#[allow(nonstandard_style)]
fn lor_from_first_vertices(vertices: &[Vertex]) -> Option<Hdf5Lor> {
    let mut in_lxe = vertices_in_scintillator(vertices);
    let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: E2, ..} = in_lxe.find(|v| v.track_id == 2)?;
    let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: E1, ..} = in_lxe.find(|v| v.track_id == 1)?;
    Some(Hdf5Lor {
        dt: t2 - t1,                   x1, y1, z1,   x2, y2, z2,
        q1: f32::NAN, q2: f32::NAN,        E1,           E2,
    })
}

#[allow(nonstandard_style)]
fn lor_from_barycentre_of_vertices(vertices: &[Vertex]) -> Option<Hdf5Lor> {
    let (a,b): (Vec<_>, Vec<_>) = vertices
        .iter()
        .filter   (|v| v.volume_id == 0 && v.parent_id <= 2)
        .partition(|v| v. track_id == 1 || v.parent_id == 1);

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

fn lor_from_hits(xyzs: &SensorMap) -> impl Fn(&[QT]) -> Option<Hdf5Lor> + '_ {
    move |hits| {
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
}

fn count_clusters(labels: &ndarray::Array1<Option<usize>>) -> usize {
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
        assert_eq!(count_clusters(&a), 0);
        let _a = array![Some(2)];
    }
}

fn lor_from_hits_dbscan(
    xyzs: &SensorMap,
    min_points: usize,
    tolerance: Length
) -> impl Fn(&[QT]) -> Option<Hdf5Lor> + '_ {
    use linfa_clustering::AppxDbscan;
    use linfa::traits::Transformer;
    move |hits| {
        let active_sensor_positions: ndarray::Array2<f32> = hits.iter()
            .flat_map(|QT { sensor_id, ..}| xyzs.get(sensor_id))
            .map(|&(x,y,z)| [mm_(x), mm_(y), mm_(z)])
            .collect::<Vec<_>>()
            .into();
        let params = AppxDbscan::params(min_points)
            .tolerance(mm_(tolerance)); // > 7mm between sipm centres
        let labels = params.transform(&active_sensor_positions).ok()?;
        let n_clusters = count_clusters(&labels);
        if n_clusters != 2 { return None }
        let mut cluster: [Vec<f32>; 2] = [vec![], vec![]];
        for (c, point) in labels.iter().zip(active_sensor_positions.outer_iter()) {
            if let Some(c) = c {
                cluster[*c].extend(point);
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
    let mut delta_E = 0.0;
    let mut rr = Length::ZERO;
    let mut xx = Length::ZERO;
    let mut yy = Length::ZERO;
    let mut zz = Length::ZERO;
    let mut tt = Time::ZERO;
    for &&Vertex { x, y, z, t, pre_KE, post_KE, .. } in vertices {
        let dE = pre_KE - post_KE;
        let (x, y, z, t) = (mm(x), mm(y), mm(z), ns(t));
        delta_E += dE;
        rr += (x*x + y*y).sqrt() * dE;
        xx += x * dE;
        yy += y * dE;
        zz += z * dE;
        // Mean time as simple compromise since want to be as close to
        // the start of event without having to much variation caused
        // by using the first detected photon (or vertex here).
        tt += t;
    };
    rr /= delta_E;
    let angle = yy.atan2(xx);
    xx  = rr * angle.cos();
    yy  = rr * angle.sin();
    tt /= vertices.len() as f32;
    Some(Barycentre { x: xx, y: yy, z: zz / delta_E, t: tt, E: delta_E })
}

#[cfg(test)]
mod test_vertex_barycentre {
    use super::*;
    use float_eq::assert_float_eq;
    use units::radian;
    use std::f32::consts::PI;

    /// Create a vertex with interesting x,y, optional pre_KE and dummy values elsewhere
    fn vertex(x: Length, y: Length, pre_ke: Option<Energyf32>) -> Vertex {
        Vertex {
            // interesting values
            x: mm_(x), y: mm_(y),
            // dummy values
            z: 23.4, t: 123.0,
            event_id: 0, parent_id: 0, track_id: 0, process_id: 0, volume_id: 0,
            moved: 0.0, deposited: 0, pre_KE: pre_ke.unwrap_or(511.0), post_KE: 0.0,
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
            .map(|(x,y)| vertex(x,y, None))
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
        let vertices = vec![vertex(mm(100.0), mm(0.0), None), vertex(mm(0.0), mm(100.0), None)];

        // The barycentre of the above points should be at 45 degrees or pi / 4.
        let angle = PI / 4.0;
        let (expected_x, expected_y) = (r * angle.cos(), r * angle.sin());

        // Create vector of vertex refs, as required by vertex_barycentre
        let vertex_refs: Vec<&Vertex> = vertices.iter().collect();

        let Barycentre { x, y, .. } = vertex_barycentre(&vertex_refs).unwrap();

        assert_float_eq!(mm_(x), expected_x, ulps <= 1);
        assert_float_eq!(mm_(y), expected_y, ulps <= 1);
    }

    #[test]
    fn energy_weights_used_correctly() {
        let energies = vec![511.0, 415.7, 350.0, 479.0, 222.5];
        let ys       = vec![353.5, 382.0, 367.3, 372.9, 377.0];
        let vertices: Vec<Vertex> = ys.iter().zip(energies.iter())
            .map(|(y, e)| vertex(mm(0.0), mm(*y), Some(*e)))
            .collect();

        // Create vector of vertex refs, as required by vertex_barycentre
        let vertex_refs: Vec<&Vertex> = vertices.iter().collect();

        let weighted_sum = ys.iter().zip(energies.iter())
            .fold(0.0, |acc, (y, w)| acc + y * w);
        let expected_e = energies.iter().sum();
        let expected_y = weighted_sum / expected_e;

        let Barycentre { x, y, E, .. } = vertex_barycentre(&vertex_refs).unwrap();

        assert_float_eq!(E, expected_e, ulps <= 1);
        assert_float_eq!(mm_(x), 0.0, abs <= 2e-5);
        assert_float_eq!(mm_(y), expected_y, ulps <=1);
    }
}

fn group_by<T>(group_by: impl FnMut(&T) -> u32, qts: impl IntoIterator<Item = T>) -> Vec<Vec<T>> {
    qts.into_iter()
        .group_by(group_by)
        .into_iter()
        .map(|(_, group)| group.collect())
        .collect()
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
        let read_data = petalo::io::hdf5::read_table::<Outer>(&file_path, "just-testing/nested", Bounds::none())?.to_vec();
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

#[cfg(test)]
mod test_dbscan {
    use linfa_clustering::{Dbscan, AppxDbscan, generate_blobs};
    use linfa::traits::Transformer;
    use ndarray::{Axis, array};
    use ndarray_rand::rand::SeedableRng;
    use rand_isaac::Isaac64Rng;

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
                arry[*c].extend(point);
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
