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


pub fn lors_from<T>(events: &[Vec<T>], mut one_lor: impl FnMut(&[T]) -> Option<Hdf5Lor>) -> Vec<Hdf5Lor> {
    events.iter()
        .flat_map(|data| one_lor(data))
        .collect()
}

pub fn filter_pde(pde: f32) -> bool {
    random::<f32>() < pde
}

pub fn gaussian_sampler(mu: f32, sigma: f32) -> Box<dyn Fn() -> f32> {
    let dist = Normal::new(mu, sigma).unwrap();
    Box::new(move || dist.sample(&mut rand::thread_rng()))
}
