/// Simulate sensor electronics using charge and time smearing
use std::path::PathBuf;

use itertools::Itertools;

use geometry::uom::ConstZero;
use geometry::units::{mm_, ns, ns_, ratio};
use geometry::units::mmps::f32::Area;
use rand::random;
use rand_distr::{Normal, Distribution};

use crate::{Length, Ratio, Time};
use crate::io::hdf5::Hdf5Lor;
use crate::io::mcreaders::{read_sensor_hits, read_sensor_map};
use crate::io::mcreaders::{MCQT, MCSensorHit, SensorMap};
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

}
