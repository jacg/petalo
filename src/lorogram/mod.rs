use ndhistogram::{ndhistogram, axis::{Uniform, Cyclic}, Histogram, HistND};

// TODO: replace with uom
type Length = f32;
type Ratio  = f32;
type Angle  = f32;

type Point = (Length, Length, Length);

/// Distinguish between true, scatter and random prompt signals
pub enum Prompt { True, Scatter, Random }

pub trait Lorogram {
    fn fill              (&mut self, p1: Point, p2: Point);
    fn value             (&    self, p1: Point, p2: Point) -> usize;
    fn interpolated_value(&    self, p1: Point, p2: Point) -> Ratio;
}

pub struct Scattergram {
    trues  : Box<dyn Lorogram>,
    scatters:Box<dyn Lorogram>,
}

impl Scattergram {

    pub fn new(make_empty_lorogram: &(dyn Fn() -> Box<dyn Lorogram>)) -> Self {
        let trues    = make_empty_lorogram();
        let scatters = make_empty_lorogram();
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
    pub fn value(&self, p1: Point, p2: Point) -> Ratio {
        let trues = self.trues.value(p1, p2);
        if trues > 0 {
            let scatters: f32 = self.scatters.value(p1, p2) as f32;
            let trues = trues as f32;
            (scatters + trues) / trues
        } else { 1.0 }
    }

    pub fn triplet(&self, p1: Point, p2: Point) -> (Ratio, f32, f32) {
        let trues = self.trues.value(p1, p2);
        if trues > 0 {
            let scatters: f32 = self.scatters.value(p1, p2) as f32;
            let trues = trues as f32;
            ((scatters + trues) / trues, trues, scatters)
        } else { (1.0, 0.0, self.scatters.value(p1, p2) as f32) }
    }
}
// --------------------------------------------------------------------------------
type Uniform1DHist = HistND<(Uniform<Length>,), usize>;

pub struct JustZ {
    histogram: Uniform1DHist,
}

impl JustZ {
    pub fn new(l: Length, nbins: usize) -> Self {
        Self { histogram: ndhistogram!(Uniform::new(nbins, -l/2.0, l/2.0); usize) }
    }
}

impl Lorogram for JustZ {

    // TODO should the points be taken by reference?

    fn fill (&mut self, p1: Point, p2: Point)          {  self.histogram.fill (&z_of_midpoint(p1, p2)); }
    fn value(&    self, p1: Point, p2: Point) -> usize { *self.histogram.value(&z_of_midpoint(p1, p2)).unwrap_or(&0) }

    fn interpolated_value(&self, _p1: Point, _p2: Point) -> f32   { todo!() }
}

fn z_of_midpoint(p1: Point, p2: Point) -> Length { (p1.2 + p2.2) / 2.0 }

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
pub struct JustR {
    histogram: Uniform1DHist,
}

impl JustR {
    pub fn new(r: Length, nbins: usize) -> Self {
        Self { histogram: ndhistogram!(Uniform::new(nbins, 0.0, r); usize) }
    }
}

impl Lorogram for JustR {
    fn fill (&mut self, p1: Point, p2: Point)          {  self.histogram.fill (&distance_from_z_axis(p1, p2));}
    fn value(&    self, p1: Point, p2: Point) -> usize { *self.histogram.value(&distance_from_z_axis(p1, p2)).unwrap_or(&0) }

    fn interpolated_value(&    self, p1: Point, p2: Point) -> Ratio {
        todo!()
    }
}

fn distance_from_z_axis(p1: Point, p2: Point) -> Length {
    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    let x1 = p1.0;
    let y1 = p1.1;
    (dx * y1 - dy * x1).abs() / (dx*dx + dy*dy).sqrt()
}
// --------------------------------------------------------------------------------
type Cyclic1DHist = HistND<(Cyclic<f32>,), usize>;

pub struct JustPhi {
    histogram: Cyclic1DHist,
}

impl JustPhi {
    pub fn new(nbins: usize) -> Self {
        Self { histogram: ndhistogram!(Cyclic::new(nbins, 0.0, std::f32::consts::PI); usize) }
    }
}

impl Lorogram for JustPhi {
    fn fill (&mut self, p1: Point, p2: Point)          {  self.histogram.fill (&phi(p1, p2)); }
    fn value(&    self, p1: Point, p2: Point) -> usize { *self.histogram.value(&phi(p1, p2)).unwrap_or(&0) }

    fn interpolated_value(&    self, p1: Point, p2: Point) -> Ratio {
        todo!()
    }
}

fn phi(p1: Point, p2: Point) -> Angle {
    let dx = p2.0 - p1.0;
    let dy = p2.1 - p1.1;
    phi_of_x_y(dx, dy)
}

fn phi_of_x_y(x: Length, y: Length) -> Angle { y.atan2(x) }
// --------------------------------------------------------------------------------
pub struct JustDeltaZ {
    histogram: Uniform1DHist,
}

impl JustDeltaZ {
    pub fn new(dz_max: Length, nbins: usize) -> Self {
        Self { histogram: ndhistogram!(Uniform::new(nbins, 0.0, dz_max); usize) }
    }
}

impl Lorogram for JustDeltaZ {
    fn fill (&mut self, p1: Point, p2: Point)          {  self.histogram.fill (&delta_z(p1, p2)); }
    fn value(&    self, p1: Point, p2: Point) -> usize { *self.histogram.value(&delta_z(p1, p2)).unwrap_or(&0) }

    fn interpolated_value(&    self, p1: Point, p2: Point) -> Ratio {
        todo!()
    }
}

fn delta_z(p1: Point, p2: Point) -> Length { (p1.2 - p2.2).abs() }
// --------------------------------------------------------------------------------
type Uniform2DHist = HistND<(Uniform<Length>, Uniform<Length>), usize>;

pub struct ZAndDeltaZ {
    histogram: Uniform2DHist
}

impl ZAndDeltaZ {
    pub fn new(l: Length, nbins_z: usize, dz_max: Length, nbins_dz: usize) -> Self {
        Self {
            histogram: ndhistogram!(
                Uniform::new(nbins_z, 0.0, dz_max),
                Uniform::new(nbins_dz, -l/2.0, l/2.0);
                usize)
        }
    }
}

impl Lorogram for ZAndDeltaZ {
    fn fill (&mut self, p1: Point, p2: Point)          {  self.histogram.fill (&(z_of_midpoint(p1, p2), delta_z(p1, p2))); }
    fn value(&    self, p1: Point, p2: Point) -> usize { *self.histogram.value(&(z_of_midpoint(p1, p2), delta_z(p1, p2))).unwrap_or(&0) }

    fn interpolated_value(&    self, p1: Point, p2: Point) -> Ratio {
        todo!()
    }
}
