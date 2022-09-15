/// Simulate sensor electronics using charge and time smearing
use std::path::PathBuf;

use itertools::Itertools;

use uom::typenum::P2;
use units::uom::ConstZero;
use units::{mm, mm_, ns, ns_, radian_, ratio, ratio_, turn};
use units::mmps::f32::Area;
use rand::random;
use rand_distr::{Normal, Distribution};

use units::{C, Quantity, Length, Ratio, Time, Angle, Velocity};
use crate::io::hdf5::Hdf5Lor;
use crate::io::mcreaders::read_sensor_hits;
use crate::io::mcreaders::{MCSensorHit, SensorMap, SensorReadout};
use crate::config::mlem::Bounds;

type DetectionEff = f32;

pub fn c_in_medium(rindex: f32) -> Velocity {
    C / ratio(rindex)
}

pub fn lors_from<T>(events: &[Vec<T>], mut one_lor: impl FnMut(&[T]) -> Option<Hdf5Lor>) -> Vec<Hdf5Lor> {
    events.iter()
        .flat_map(|data| one_lor(data))
        .collect()
}

pub struct LorBatch {
    pub lors: Vec<Hdf5Lor>,
    pub n_events_processed: usize,
}

impl LorBatch {
    pub fn new(lors: Vec<Hdf5Lor>, n_events_processed: usize) -> Self {
        Self { lors, n_events_processed }
    }
}

pub fn random_detection_pde<T>(pde: DetectionEff) -> impl Fn(&T) -> bool {
    move |_| random::<DetectionEff>() < pde
}

pub fn gaussian_sampler(mu: f32, sigma: f32) -> impl Fn() -> f32 {
    let dist = Normal::new(mu, sigma).unwrap();
    move || dist.sample(&mut rand::thread_rng())
}

// Could this be a FROM? How to use the Group here to avoid too many collects?
fn qt_min_time(sensor_id: u32, hits: Vec<MCSensorHit>, xyzs: &SensorMap) -> Option<SensorReadout> {
    let first_time = hits.iter().fold(f32::INFINITY, |min_val, val| val.time.min(min_val));
    let &(x, y, z) = xyzs.get(&sensor_id)?;
    Some(SensorReadout { sensor_id, x, y, z, q: hits.len() as u32, t: ns(first_time) })
}

/// combine the individual hits into information for each sensor.
/// Currently uses the minimum as the time for each sensor which
/// might prove too inflexible. Should be generalised.
pub fn combine_sensor_hits(hits: Vec<MCSensorHit>, xyzs: &SensorMap) -> Vec<SensorReadout> {
    hits.into_iter()
        .sorted_by_key(|&MCSensorHit { sensor_id, ..}| sensor_id)
        .group_by(|h| h.sensor_id)
        .into_iter()
        .flat_map(|(s, grp)| qt_min_time(s, grp.collect(), &xyzs))
        .collect()
}

// This seems like something that would be in a crate somewhere.
fn dot((x1,y1,z1): (Length, Length, Length), (x2,y2,z2): (Length, Length, Length)) -> Area {
    x1*x2 + y1*y2 + z1*z2
}

// Dot product separation for now, need to implement other algorithm(s)
fn dot_product_clusters(hits: &[SensorReadout], threshold: f32, min_sensors: usize) -> Option<(Vec<SensorReadout>, Vec<SensorReadout>)> {
    let max_pos = hits.iter().max_by_key(|e| e.q).map(|&SensorReadout { x, y, z, .. }| (x, y, z))?;
    let (c1, c2): (Vec<SensorReadout>, Vec<SensorReadout>) = hits.iter()
            .filter(|&SensorReadout { q, .. }| *q as f32 > threshold)
            .partition(|&SensorReadout { x, y, z, .. }| dot(max_pos, (*x, *y, *z)) > Area::ZERO);
    if c1.len() >= min_sensors && c2.len() >= min_sensors { Some((c1, c2)) } else { None }
}

// Is there a more generic way to avoid rewriting this code over and over again?
fn get_sensor_charge(hits: &[SensorReadout]) -> Vec<Ratio> {
    hits.iter()
        .map(|&SensorReadout { q, .. }| ratio(q as f32))
        .collect()
}

fn get_sensor_z(hits: &[SensorReadout]) -> Vec<Length> {
    hits.iter()
        .map(|&SensorReadout { z, .. }| z)
        .collect()
}

/// Reconstruct LORs based on the barycentre of the 2 selected clusters (dot prod).
fn lor_from_sensor_positions(
    hits        : &[SensorReadout],
    threshold   : f32             ,
    min_sensors : usize           ,
    min_charge  : f32             ,
    doi_function: impl Fn(Barycentre, &[SensorReadout]) -> Option<Barycentre>
) -> Option<Hdf5Lor> {
    if hits.is_empty() { return None }
    let (c1, c2) = dot_product_clusters(hits, threshold, min_sensors)?;
    // TODO time just taken as minimum of minima. Averaging used elsewhere.\
    let b1 = sipm_charge_barycentre(&c1);
    let b1 = doi_function(b1, &c1)?;
    let b2 = sipm_charge_barycentre(&c2);
    let b2 = doi_function(b2, &c2)?;
    if b1.q.min(b2.q) < ratio(min_charge) { return None }
    Some(Hdf5Lor {
        dt: ns_(b2.t - b1.t),
        x1: mm_(b1.x), y1: mm_(b1.y), z1: mm_(b1.z),
        x2: mm_(b2.x), y2: mm_(b2.y), z2: mm_(b2.z),
        q1: ratio_(b1.q), q2: ratio_(b2.q), E1: f32::NAN, E2: f32::NAN,
    })
}

#[derive(Copy, Clone)]
struct Barycentre {
    x: Length,
    y: Length,
    z: Length,
    t: Time,
    q: Ratio,
}

// TODO protect against curvature reducing R for big clusters?
fn sipm_charge_barycentre(hits: &[SensorReadout]) -> Barycentre {
    let mut qs = Ratio::ZERO;
    let mut xx = Length::ZERO;
    let mut yy = Length::ZERO;
    let mut zz = Length::ZERO;
    let mut tt = Time::ZERO;
    for &SensorReadout{ x, y, z, q, t, .. } in hits {
        let q = ratio(q as f32);
        qs += q;
        xx += x * q;
        yy += y * q;
        zz += z * q;
        if t < tt { tt = t };
    }
    Barycentre { x: xx / qs, y: yy / qs, z: zz / qs, t: tt, q: qs }
}

pub fn lor_reconstruction<'a>(
    xyzs       : &'a SensorMap,
    pde        : f32,
    time_mu    : f32,
    time_sigma : f32,
    threshold  : f32,
    min_sensors: usize,
    min_charge : f32,
) -> Box<dyn Fn(&PathBuf) -> hdf5::Result<LorBatch> + 'a> {
    let time_smear = gaussian_sampler(time_mu, time_sigma);
    // TODO Values should be configurable, need to generalise.
    let doi_func = calculate_interaction_position(DOI::ZRMS, ratio(-1.2906), mm(384.4280), mm(352.0), mm(382.0));
    Box::new(move |filename: &PathBuf| -> hdf5::Result<LorBatch> {
        let sensor_hits = read_sensor_hits(&filename, Bounds::none())?;
        let detected_hits: Vec<MCSensorHit> =
            sensor_hits.iter()
                       .filter(random_detection_pde(pde))
                       .map(|hit| MCSensorHit { time: hit.time + time_smear(), ..*hit })
                       .collect();
        let events: Vec<Vec<SensorReadout>> =
            detected_hits.into_iter()
                         .group_by(|h| h.event_id)
                         .into_iter()
                         .map(|(_, grp)| combine_sensor_hits(grp.collect(), &xyzs))
                         .collect();
        Ok(LorBatch::new(lors_from(&events, |evs| lor_from_sensor_positions(evs, threshold, min_sensors, min_charge, &doi_func)), events.len()))
    })
}

enum DOI { ZRMS, PRMS, CRMS }

// Function to correct for DOI
// Takes into account the physical range for the setup.
// TODO Only implemented for zrms, need phi and combined
fn calculate_interaction_position(
    ftype: DOI   ,
    m    : Ratio ,
    c    : Length,
    rmin : Length,
    rmax : Length
) -> Box<dyn Fn(Barycentre, &[SensorReadout]) -> Option<Barycentre>> {
    match ftype {
        DOI::ZRMS =>
            Box::new(move |barycentre: Barycentre, hits: &[SensorReadout]| -> Option<Barycentre> {
                let zs = get_sensor_z(hits);
                let weights = get_sensor_charge(hits);
                let rms = weighted_std(&zs, &weights, Some(barycentre.z))?;
                let mut r = m * rms + c;// Bias correction?
                if r < rmin { r = rmin; }
                else if r > rmax { r = rmax; }
                let phi = barycentre.y.atan2(barycentre.x);
                let t = correct_interaction_time(&barycentre, rmax - r);
                Some(Barycentre {x: r * phi.cos(), y: r * phi.sin(), t, ..barycentre})
            }),
        DOI::PRMS =>
            Box::new(move |barycentre: Barycentre, hits: &[SensorReadout]| -> Option<Barycentre> {
                let phis = azimuthal_angle(hits)?;
                let weights = get_sensor_charge(hits);
                let rms = weighted_std(&phis, &weights, None)?;
                // The parameterisation here isn't strictly units safe as use a std of angle to get a length.
                let mut r = mm(ratio_(m) * radian_(rms)) + c;
                if r < rmin { r = rmin; }
                else if r > rmax { r = rmax; }
                let phi = barycentre.y.atan2(barycentre.x);
                let t = correct_interaction_time(&barycentre, rmax - r);
                Some(Barycentre {x: r * phi.cos(), y: r * phi.sin(), t, ..barycentre})
            }),
        DOI::CRMS => 
            Box::new(move |barycentre: Barycentre, hits: &[SensorReadout]| -> Option<Barycentre> {
                let zs = get_sensor_z(hits);
                let phis = azimuthal_angle(hits)?;
                let weights = get_sensor_charge(hits);
                let zrms: Length = weighted_std(&zs, &weights, Some(barycentre.z))?;
                let prms: Angle = weighted_std(&phis, &weights, None)?;
                let mut r = m * (zrms.powi(P2::new()) + (rmax * prms).powi(P2::new())).sqrt() + c;
                if r < rmin { r = rmin; }
                else if r > rmax { r = rmax; }
                let phi = barycentre.y.atan2(barycentre.x);
                let t = correct_interaction_time(&barycentre, rmax - r);
                Some(Barycentre {x: r * phi.cos(), y: r * phi.sin(), t, ..barycentre})
            }),
    }
}

fn correct_interaction_time(b: &Barycentre, rdiff: Length) -> Time {
    b.t - rdiff / c_in_medium(1.69)
}

// Calculate the azimuthal angle for each hit and
// adjust so continuous (make option?)
fn azimuthal_angle(hits: &[SensorReadout]) -> Option<Vec<Angle>> {
    if hits.is_empty() { return None }
    let mut phis: Vec<Angle> = hits.iter()
        .map(|&SensorReadout { x, y, .. }| y.atan2(x))
        .collect();
    
    if phis.iter().any(|&phi| phi < turn(-0.25)) && phis.iter().any(|&phi| phi > turn(0.0)) {
        // There is a discontinuity in this case (sensors in second and third quadrants) so mean or std calculations would be wrong.
        for phi in phis.iter_mut().filter(|p| **p < turn(-0.25)) {
            *phi += turn(1.0);
        }
    }
    Some(phis)
}

// Weighted mean and variance.
pub fn weighted_mean<D1, D2, U>(data: &[Quantity<D1, U, f32>], weights: &[Quantity<D2, U, f32>]) -> Option<Quantity<D1, U, f32>>
where
    D1: uom::si::Dimension  + ?Sized,
    D2: uom::si::Dimension  + ?Sized,
    U : uom::si::Units<f32> + ?Sized,
{
    if data.is_empty() { return None }
    assert_eq!(data.len(), weights.len());
    assert!(weights.iter().all(|&w| w > Quantity::<D2, U, f32>::ZERO));

    let wsum: f32 = weights.iter().cloned().map(|w| w.value).sum();
    let weighted_sum: f32 = data.iter().cloned()
        .zip(weights.iter().cloned())
        .map(|(d, w)| (d.value, w.value))
        .map(|(d, w)| w * d)
        .sum();
    Some(
        Quantity::<D1, U, f32> {
            dimension: uom::lib::marker::PhantomData,
            units: uom::lib::marker::PhantomData,
            value: weighted_sum / wsum,
        }
    )
}

pub fn weighted_std<D1, D2, U>(data: &[Quantity<D1, U, f32>], weights: &[Quantity<D2, U, f32>], mean: Option<Quantity<D1, U, f32>>) -> Option<Quantity<D1, U, f32>>
where
    D1: uom::si::Dimension  + ?Sized,
    D2: uom::si::Dimension  + ?Sized,
    U : uom::si::Units<f32> + ?Sized,
{
    if data.is_empty() { return None }
    assert_eq!(data.len(), weights.len());
    assert!(weights.iter().all(|&w| w > Quantity::<D2, U, f32>::ZERO));

    let uom_mu = match mean {
        Some(mval) => mval,
        None => weighted_mean(data, weights)?
    };
    let mu = uom_mu.value;

    let adj_wsum: f32 = weights.iter().cloned().map(|w| w.value).sum::<f32>() - 1.0;
    let var: f32 = data.iter().cloned()
        .zip(weights.iter().cloned())
        .map(|(d, w)| (d.value, w.value))
        .map(|(d, w)| w * (d - mu).powi(2))
        .sum();
    Some(
        Quantity::<D1, U, f32> {
            dimension: uom::lib::marker::PhantomData,
            units: uom::lib::marker::PhantomData,
            value: (var / adj_wsum).sqrt(),
        }
    )
}


// ----- Tests -----

#[cfg(test)]
mod test_electronics {
    use super::*;
    use rstest::{fixture, rstest};
    use float_eq::assert_float_eq;

    #[fixture]
    fn mcsensorhit_vector() -> ([u32; 5], SensorMap, Vec<MCSensorHit>) {
        let dummy_times: [f32; 5] = [0.4, 0.25, 0.11, 0.50, 0.47];
        let sensor_ids: [u32; 5]  = [133, 133, 133, 135, 136];
        let positions = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                         (0.0, 0.0, 0.0), (0.0, 0.0, 1.0),
                         (0.0, 1.0, 0.0)];
        let smap = sensor_ids.iter().cloned()
            .zip(positions.iter())
            .map(|(id, &(x, y, z))| (id, (mm(x), mm(y), mm(z))))
            .collect();
        let event_id: u32 = 1000;
        let dummy_hits = dummy_times
            .iter()
            .zip(sensor_ids.iter())
            .map(|(&t, &s)| MCSensorHit { event_id, sensor_id: s, time: t})
            .collect();
        (sensor_ids, smap, dummy_hits)
    }

    // Some dummy sensors to test an easy separation.
    #[fixture]
    fn qt_vector() -> ([u32; 10], Vec<SensorReadout>) {
        let sensor_ids: [u32; 10] = [0, 1, 2, 3, 4, 5, 20, 21, 22, 23];
        let positions = [(0.0,  0.0), ( 0.0, 1.0),
                         (1.0,  1.0), ( 1.0, 0.0),
                         (1.0, -1.0), (-1.0, 0.0),
                         (0.0,  0.0), ( 0.0, 1.0),
                         (1.0,  1.0), ( 1.0, 0.0)];
        let (p1, p2): (f32, f32) = (-20.0, 20.0);
        let qs: [u32; 10] = [5, 3, 2, 3, 1, 2, 7, 2, 3, 2];
        let ts: [f32; 10] = [0.11; 10];
        let qts = sensor_ids.iter()
            .zip(positions.iter())
            .zip(qs.iter())
            .zip(ts.iter())
            .map(|(((&id, &pos), &q), &t)| {
                let p = if id < 10 { p1 } else {p2};
                SensorReadout { sensor_id: id, x: mm(p), y: mm(p + pos.0), z: mm(p + pos.1), q, t: ns(t) }
            })
            .collect();
        (qs, qts)
    }

    #[test]
    fn test_random_detection_pde() {
        let pde: f32 = 0.45;

        let n = 1000;
        let selected: Vec<i32> = (0..n).filter(random_detection_pde(pde)).collect();
        let samp_3sig = 3.0 * (0.45 * (1.0 - 0.45) / n as f32).sqrt();
        assert_float_eq!(selected.len() as f32 / n as f32, pde, abs <= samp_3sig)
    }

    #[rstest]
    fn test_sensor_min(mcsensorhit_vector: ([u32; 5], SensorMap, Vec<MCSensorHit>)) {
        let (sensor_ids, smap, hits) = mcsensorhit_vector;

        let qt = qt_min_time(sensor_ids[0], hits[0..3].to_vec(), &smap).unwrap();
        assert_eq!(qt.q, 3);
        assert_float_eq!(ns_(qt.t), 0.11, ulps <=1)
    }

    #[rstest]
    fn test_combine_sensors(mcsensorhit_vector: ([u32; 5], SensorMap, Vec<MCSensorHit>)) {
        let (_, smap, hits) = mcsensorhit_vector;

        let qts = combine_sensor_hits(hits, &smap);
        assert_eq!(qts.len(), 3);
        let (qs, ts): (Vec<u32>, Vec<f32>) = qts.iter().map(|&SensorReadout { q, t, .. }| (q, ns_(t))).unzip();
        assert_eq!(qs, [3, 1, 1]);
        // Probably want to check all the ts too but getting compilation errors. TODO
        assert_float_eq!(ts[0], 0.11, ulps <= 1)
    }

    #[test]
    fn test_combine_sensors_unordered() {
        let sensor_ids = [132, 145, 129, 129, 132, 145];
        let positions: [(f32, f32, f32); 6] = [(0.0, 0.0, 1.0), (1.0, 1.0, 1.0),
                                               (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                               (0.0, 0.0, 1.0), (1.0, 1.0, 1.0)];
        let smap = sensor_ids.iter().cloned()
            .zip(positions.iter())
            .map(|(id, &(x, y, z))| (id, (mm(x), mm(y), mm(z))))
            .collect();
        let test_hits: Vec<MCSensorHit> = sensor_ids
            .iter()
            .map(|&sensor_id| MCSensorHit { event_id: 1, sensor_id, time: 0.11})
            .collect();
        let qts = combine_sensor_hits(test_hits, &smap);
        assert_eq!(qts.len(), 3);
        let ids: Vec<u32> = qts
            .iter()
            .map(|SensorReadout { sensor_id, .. }| *sensor_id)
            .collect();
        assert_eq!(ids, [129, 132, 145]);
    }

    #[rstest]
    fn test_dotproduct_separation(qt_vector: ([u32; 10], Vec<SensorReadout>)) {
        let (qs, sensors) = qt_vector;
        let clustering = dot_product_clusters(&sensors, 1.0, 2);
        assert!(clustering.is_some());
        let (c1, c2) = clustering.unwrap();
        assert_eq!(c1.len(), 4);
        assert_eq!(c2.len(), 5);
        assert_eq!(c1.iter().map(|&SensorReadout { q, .. }| q).sum::<u32>(), qs[6..10].iter().sum());
        assert_eq!(c2.iter().map(|&SensorReadout { q, .. }| q).sum::<u32>(), qs[0..6].iter().sum::<u32>() - 1);
    }

    #[rstest]
    fn test_lor_from_sensor_position(qt_vector: ([u32; 10], Vec<SensorReadout>)) {
        let (_, sensors) = qt_vector;
        let lor = lor_from_sensor_positions(&sensors, 1.0, 2, 2.0,
            |b: Barycentre, _h: &[SensorReadout]| Some(b));
        assert!(lor.is_some());
    }

    use units::assert_uom_eq;
    use uom::si::length::millimeter;
    use uom::si::angle::radian;
    #[rstest]
    fn test_length_weighted_mean(qt_vector: ([u32; 10], Vec<SensorReadout>)) {
        let (qs, sensors) = qt_vector;
        let mut xs = Vec::with_capacity(qs.len());
        let mut qtest = Vec::with_capacity(qs.len());
        let mut xq_sum = Length::ZERO * Ratio::ZERO;
        let mut wsum = Ratio::ZERO;
        for (&SensorReadout { x, .. }, &q) in sensors.iter().zip(qs.iter()) {
            xs.push(x);
            let qrat = ratio(q as f32);
            qtest.push(qrat);
            xq_sum += x * qrat;
            wsum += qrat;
        }
        let ave = weighted_mean(&xs, &qtest).unwrap();
        assert_uom_eq!(millimeter, ave, xq_sum / wsum, ulps <= 1)
    }

    #[rstest]
    fn test_length_weighted_std(qt_vector: ([u32; 10], Vec<SensorReadout>)) {
        let (qs, sensors) = qt_vector;
        let mut xs = Vec::with_capacity(qs.len());
        let mut qtest = Vec::with_capacity(qs.len());
        let mut wsum = Ratio::ZERO;
        let mut ave = Length::ZERO;
        let mut std = Area::ZERO * Ratio::ZERO;
        for (&SensorReadout { x, .. }, &q) in sensors.iter().zip(qs.iter()) {
            xs.push(x);
            let qrat = ratio(q as f32);
            qtest.push(qrat);
            wsum += qrat;
            let mean_old = ave;
            ave = mean_old + (qrat / wsum) * (x - mean_old);
            std += qrat * (x - mean_old) * (x - ave);
        }
        let adj_wsum = wsum - ratio(1.0);
        let std_calc = weighted_std(&xs, &qtest, None).unwrap();
        let std_expt = (std / adj_wsum).sqrt();
        assert_uom_eq!(millimeter, std_calc, std_expt, ulps <= 1);
        let std_calc = weighted_std(&xs, &qtest, Some(ave)).unwrap();
        assert_uom_eq!(millimeter, std_calc, std_expt, ulps <= 1);
    }

    use units::TWOPI;
    #[rstest]
    fn test_azimuthal_angle(qt_vector: ([u32; 10], Vec<SensorReadout>)) {
        let (_, hits) = qt_vector;
        let phis1 = azimuthal_angle(&hits[..6]);
        let phis_all = azimuthal_angle(&hits);
        assert!(phis1.is_some() && phis_all.is_some());
        assert!(phis1.clone().unwrap().iter().all(|phi| *phi < turn(0.0)));
        assert!(phis_all.clone().unwrap().iter().all(|phi| *phi > turn(0.0)));
        for (&p1, &p2) in phis1.unwrap().iter().zip(phis_all.unwrap().iter()) {
            assert_uom_eq!(radian, p1, p2 - TWOPI, ulps <= 1);
        }
    }

    #[rstest]
    fn test_angle_weighted_std(qt_vector: ([u32; 10], Vec<SensorReadout>)) {
        let (qs, hits) = qt_vector;
        let mut qrs = Vec::with_capacity(6);
        let mut phis = Vec::with_capacity(6);
        let mut wsum = Ratio::ZERO;
        let mut ave = Angle::ZERO;
        let mut std = Angle::ZERO;
        for (&SensorReadout { x, y, .. }, &q) in hits.iter().take(6).zip(qs.iter()) {
            let phi: Angle = y.atan2(x);
            phis.push(phi);
            let qrat = ratio(q as f32);
            qrs.push(qrat);
            wsum += qrat;
            let mean_old = ave;
            ave = mean_old + Angle::from((qrat / wsum) * (phi - mean_old));
            std += Angle::from(qrat * (phi - mean_old) * (phi - ave));
        }
        let phi_std = weighted_std(&phis, &qrs, None);
        assert!(phi_std.is_some());
        let adj_wsum = wsum - ratio(1.0);
        let std_expt: Angle = (std / adj_wsum).sqrt().into();
        // Why less accurate?
        assert_uom_eq!(radian, phi_std.unwrap(), std_expt, ulps <= 20);
        let phi_std = weighted_std(&phis, &qrs, Some(ave)).unwrap();
        assert_uom_eq!(radian, phi_std, std_expt, ulps <= 20);
    }
}
