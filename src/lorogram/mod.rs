pub mod axis;

use ndhistogram::{ndhistogram, axis::Uniform, Histogram, HistND};
use std::f32::consts::TAU;
use axis::Cyclic;

type Point = (f32, f32, f32);

/// Distinguish between true, scatter and random prompt signals
pub enum Prompt { True, Scatter, Random }

pub trait Lorogram {
    fn fill              (&mut self, p1: Point, p2: Point);
    fn value             (&    self, p1: Point, p2: Point) -> usize;
    fn interpolated_value(&    self, p1: Point, p2: Point) -> f32;
}

pub struct Scattergram<T: Lorogram> {
    trues: T,
    scatters: T,
}

impl<T: Lorogram + Clone> Scattergram<T> {
    pub fn new(lorogram: T) -> Self {
        let trues = lorogram;
        let scatters = trues.clone();
        Self { trues, scatters }
    }

    pub fn fill(&mut self, kind: Prompt, p1: Point, p2: Point) {
        match kind {
            Prompt::True    => self.trues.   fill(p1, p2),
            Prompt::Scatter => self.scatters.fill(p1, p2),
            Prompt::Random  => panic!("Not expecting any random events yet."),
        }
    } 

    /// Multiplicative contribution of scatters to trues, in nearby LORs.
    ///
    /// `(scatters + trues) / trues`
    pub fn value(&self, p1: Point, p2: Point) -> f32 {
        let trues = self.trues.value(p1, p2);
        if trues > 0 {
            let scatters: f32 = self.scatters.value(p1, p2) as f32;
            let trues = trues as f32;
            (scatters + trues) / trues
        } else { 1.0 }
    }

}

#[derive(Clone)]
pub struct JustZ {
    histogram: JustZHist,
}

type JustZHist = HistND<(Uniform<f32>,), usize>;

impl JustZ {
    pub fn new(l: f32, nbins: usize) -> Self {
        Self { histogram: ndhistogram!(Uniform::new(nbins, -l/2.0, l/2.0); usize) }
    }
}

impl Lorogram for JustZ {

    // TODO should the points be taken by reference?

    fn fill(&mut self, p1: Point, p2: Point) {
        let z = (p1.2 + p2.2) / 2.0;
        self.histogram.fill(&z);
    }

    fn value(&self, p1: Point, p2: Point) -> usize {
        let z = (p1.2 + p2.2) / 2.0;
        *self.histogram.value(&z).unwrap_or(&0)
    }

    fn interpolated_value(&self, p1: Point, p2: Point) -> f32   { todo!() }
}

#[cfg(test)]
mod test_just_z {
    use super::*;

    #[test]
    fn retrieve() {
        let mut lg = JustZ::new(1000.0, 10);
        lg.fill         ((0.0, 0.0, 111.0), (0.0, 0.0, 555.0));
        let n = lg.value((1.0, 2.0, 222.0), (9.0, 8.0, 444.0));
        assert_eq!(n, 1);
    }
}

