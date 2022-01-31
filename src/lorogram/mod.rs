pub mod axis;

use ndhistogram::{ndhistogram, axis::Uniform, Histogram, HistND};

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
// --------------------------------------------------------------------------------
type JustZHist = HistND<(Uniform<f32>,), usize>;

#[derive(Clone)]
pub struct JustZ {
    histogram: JustZHist,
}

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
// --------------------------------------------------------------------------------
type ZandTanThetaHist = HistND<(Uniform<f32>, Uniform<f32>), usize>;

#[derive(Clone)]
struct ZandTanTheta {
    histogram: ZandTanThetaHist,
}

// impl Lorogram for ZandTanTheta {
//     fn fill(&mut self, p1: Point, p2: Point) {
//         let ((x1,y1,z1), (x2,y2,z2)) = (p1,p2);
//         let z = (z1 + z2) / 2.0;
//     }
// }

// #[derive(Clone)]
// pub struct ZandTanTheta {

// }

// type Aa = HistND<(Uniform<f32>,   // lower z
//                   Uniform<f32>,   // higher z
//                   Cyclic <f32>,   // angle of lower z
//                   Cyclic <f32>,), // angle of higher z
//                   usize>;


// pub struct Aaa {
//     l: f32,
//     r: f32,
//     histogram: Aa,
// }

// impl Aaa {
//     fn new(r: f32, l: f32, n_axial_bins: usize, n_angluar_bins: usize) -> Self {
//         Self {
//             l, r,
//             histogram: ndhistogram!(
//                 Uniform::new(n_axial_bins  , -l/2.0, l/2.0), // lower z
//                 Uniform::new(n_axial_bins  , -l/2.0, l/2.0), // higher z
//                 Cyclic ::new(n_angluar_bins,    0.0, TAU  ), // phi  : around z, anticlockwise from x
//                 Cyclic ::new(n_angluar_bins,    0.0, TAU  ); // theta: around r, anticlockwise from perpendicular to z
//                 usize
//             )
//         }
//     }
// }

// impl Lorogram for Aaa {
//     fn fill(&mut self, p1: Point, p2: Point) {
//         // Ensure that p1 has the lower z-coordinate
//         let (p1, p2) = if p1.2 < p2.2 { (p1, p2) } else { (p2, p1) };
//     }
//     fn value             (&self, p1: Point, p2: Point) -> usize { todo!() }
//     fn interpolated_value(&self, p1: Point, p2: Point) -> f32   { todo!() }
// }
