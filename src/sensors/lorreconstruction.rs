/// Simulate sensor electronics using charge and time smearing
use std::path::PathBuf;

use geometry::Quantity;
use itertools::Itertools;

use geometry::uom::ConstZero;
use geometry::units::{mm_, ns, ns_, ratio};
use geometry::units::mmps::f32::Area;
use rand::random;
use rand_distr::{Normal, Distribution};

use crate::{C, Length, Ratio, Time, Angle, Velocity};
use crate::io::hdf5::Hdf5Lor;
use crate::io::mcreaders::read_sensor_hits;
use crate::io::mcreaders::{MCSensorHit, SensorMap, SensorReadout};
use crate::config::mlem::Bounds;

type DetectionEff = f32;


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

    use geometry::assert_uom_eq;
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
