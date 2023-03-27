#[allow(nonstandard_style)]
pub (super) fn lor_from_first_vertices(vertices: &[Vertex]) -> Option<Hdf5Lor> {
    let mut in_lxe = vertices_in_scintillator(vertices);
    let Vertex{x:x2, y:y2, z:z2, t:t2, pre_KE: E2, ..} = in_lxe.find(|v| v.track_id == 2)?;
    let Vertex{x:x1, y:y1, z:z1, t:t1, pre_KE: E1, ..} = in_lxe.find(|v| v.track_id == 1)?;
    Some(Hdf5Lor {
        dt: t2 - t1,                   x1, y1, z1,   x2, y2, z2,
        q1: f32::NAN, q2: f32::NAN,        E1,           E2,
    })
}

#[allow(nonstandard_style)]
pub (super) fn lor_from_barycentre_of_vertices(vertices: &[Vertex]) -> Option<Hdf5Lor> {
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

pub (super) fn lor_from_hits(xyzs: &SensorMap) -> impl Fn(&[QT]) -> Option<Hdf5Lor> + '_ {
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

pub (super) fn lor_from_hits_dbscan(
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

// ----- Imports -----------------------------------------------------------------------------------------
use petalo::{
    Point,
    io::hdf5::{Hdf5Lor, mc::Vertex},
};
use crate::{QT, SensorMap, vertices_in_scintillator};
use units::{
    Area, Length, Ratio, Time,
    todo::Energyf32,
    mm, mm_, ns, ns_, ratio,
    uom::ConstZero,
};

// ----- TESTS ------------------------------------------------------------------------------------------
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
