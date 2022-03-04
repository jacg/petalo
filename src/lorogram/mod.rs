use ndhistogram::{ndhistogram, axis::{Axis, Uniform, Cyclic}, Histogram, HistND};
use crate::weights::LOR;
use crate::types::Point;
use std::f32::consts::PI;

// TODO: replace with uom
type Length = f32;
type Ratio  = f32;
type Angle  = f32;

/// Distinguish between true, scatter and random prompt signals
pub enum Prompt { True, Scatter, Random }

pub struct Scattergram {
    trues  : Box<dyn LorogramBis>,
    scatters:Box<dyn LorogramBis>,
}

impl Scattergram {

    pub fn new(make_empty_lorogram: &(dyn Fn() -> Box<dyn LorogramBis>)) -> Self {
        let trues    = make_empty_lorogram();
        let scatters = make_empty_lorogram();
        Self { trues, scatters }
    }

    pub fn fill(&mut self, kind: Prompt, lor: &LOR) {
        match kind {
            Prompt::True    => self.trues.   put(lor),
            Prompt::Scatter => self.scatters.put(lor),
            Prompt::Random  => panic!("Not expecting any random events yet."),
        }
    } 

    /// Multiplicative contribution of scatters to trues, in nearby LORs.
    ///
    /// `(scatters + trues) / trues`
    pub fn value(&self, lor: &LOR) -> Ratio {
        let trues = self.trues.get(lor);
        if trues > 0 {
            let scatters: f32 = self.scatters.get(lor) as f32;
            let trues = trues as f32;
            (scatters + trues) / trues
        } else { 1.0 }
    }

    pub fn triplet(&self, lor: &LOR) -> (Ratio, f32, f32) {
        let trues = self.trues.get(lor);
        if trues > 0 {
            let scatters: f32 = self.scatters.get(lor) as f32;
            let trues = trues as f32;
            ((scatters + trues) / trues, trues, scatters)
        } else { (1.0, 0.0, self.scatters.get(lor) as f32) }
    }
}
// --------------------------------------------------------------------------------
pub struct MappedAxis<T,A>
where
    A: Axis,
{
    axis: A,
    map: Box<dyn Fn(&T) -> A::Coordinate>,
}

impl<T,A> Axis for MappedAxis<T,A>
where
    A: Axis,
{
    type Coordinate = T;

    type BinInterval = A::BinInterval;

    fn index(&self, coordinate: &Self::Coordinate) -> Option<usize> {
        self.axis.index(&(self.map)(coordinate))
    }

    fn num_bins(&self) -> usize {
        self.axis.num_bins()
    }

    fn bin(&self, index: usize) -> Option<Self::BinInterval> {
        self.axis.bin(index)
    }
}
// --------------------------------------------------------------------------------
pub type LorAxU = MappedAxis<LOR, Uniform<Length>>;
pub type LorAxC = MappedAxis<LOR, Cyclic <Length>>;
pub type LorogramU   = HistND<(LorAxU,               ), usize>;
pub type LorogramC   = HistND<(LorAxC,               ), usize>;
pub type LorogramUU  = HistND<(LorAxU, LorAxU        ), usize>;
pub type LorogramUC  = HistND<(LorAxU, LorAxC        ), usize>;
pub type LorogramUUU = HistND<(LorAxU, LorAxU, LorAxU), usize>;
pub type LorogramUUC = HistND<(LorAxU, LorAxU, LorAxC), usize>;


fn z_of_midpoint(LOR {p1, p2, ..}: &LOR) -> Length { (p1.z + p2.z) / 2.0 }

fn delta_z(LOR{p1, p2, ..}: &LOR) -> Length { (p1.z - p2.z).abs() }

fn distance_from_z_axis(LOR{ p1, p2, .. }: &LOR) -> Length {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let x1 = p1.x;
    let y1 = p1.y;
    (dx * y1 - dy * x1).abs() / (dx*dx + dy*dy).sqrt()
}

fn phi(LOR{ p1, p2, .. }: &LOR) -> Angle {
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    phi_of_x_y(dx, dy)
}

fn phi_of_x_y(x: Length, y: Length) -> Angle { y.atan2(x) }

pub fn axis_z(nbins: usize, min: Length, max: Length) -> LorAxU {
    LorAxU {
        axis: Uniform::new(nbins, min, max),
        map: Box::new(z_of_midpoint),
    }
}

pub fn axis_delta_z(nbins: usize, max: Length) -> LorAxU {
    LorAxU {
        axis: Uniform::new(nbins, 0.0, max),
        map: Box::new(delta_z),
    }
}

pub fn axis_r(nbins: usize, max: Length) -> LorAxU {
    LorAxU {
        axis: Uniform::new(nbins, 0.0, max),
        map: Box::new(distance_from_z_axis),
    }
}

pub fn axis_phi(nbins: usize) -> LorAxC {
    LorAxC {
        axis: Cyclic::new(nbins, 0.0, PI),
        map: Box::new(phi),
    }
}

#[cfg(test)]
mod test_mapped_axes {
    use super::*;

    #[test]
    fn uniform() {
        let nbins = 10;
        let axis = axis_phi(nbins);
        assert_eq!(axis.num_bins(), nbins);
        let mut h = ndhistogram!(axis; usize);
        let x = 150.0;
        let y = 234.5;
        let (dummy1, dummy2, dummy3, dummy4) = (111.1, 222.2, 333.3, 444.4);
        let (a, b) = (30.0, 40.0); // scaling factors
        h.fill         (&LOR::from(((a*x, a*y, dummy1), (-a*x, -a*y, dummy2))));
        let n = h.value(&LOR::from(((b*x, b*y, dummy3), (-b*x, -b*y, dummy4)))).unwrap_or(&0);
        assert_eq!(n, &1);
    }

    #[test]
    fn two_dimensions() {
        let nbins_z = 10;
        let nbins_dz = 10;
        let l = 1000.0;
        let max_dz = l;
        let mut h = ndhistogram!(
            axis_z      (nbins_z , -l/2.0, l/2.0),
            axis_delta_z(nbins_dz, max_dz);
            usize
        );
        let (z, delta) = (123.4, 543.2);
        // Irrelevant values
        let (i1, i2, i3, i4, i5, i6, i7, i8) = (10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0);

        let l1 = LOR::from(((i1, i2, z-delta), (i3, i4, z+delta)));
        let l2 = LOR::from(((i5, i6, z+delta), (i7, i8, z-delta)));
        h.put        (&l1);
        let n = h.get(&l2);

        assert_eq!(n, 1);
    }

}
// --------------------------------------------------------------------------------
pub trait LorogramBis {
    fn put(&mut self, lor: &LOR);
    fn get(&    self, lor: &LOR) -> usize;
}

impl<X> LorogramBis for ndhistogram::Hist1D<X, usize>
where
    X: Axis<Coordinate = LOR>,
{
    fn put(&mut self, lor: &LOR)          {  self.fill (lor) }
    fn get(&    self, lor: &LOR) -> usize { *self.value(lor).unwrap_or(&0) }
}

impl<X, Y> LorogramBis for ndhistogram::Hist2D<X, Y, usize>
where
    X: Axis<Coordinate = LOR>,
    Y: Axis<Coordinate = LOR>,
{
    fn put(&mut self, lor: &LOR)          {  self.fill (&(*lor, *lor)) }
    fn get(&    self, lor: &LOR) -> usize { *self.value(&(*lor, *lor)).unwrap_or(&0) }
}
// --------------------------------------------------------------------------------
impl From<((f32, f32, f32), (f32, f32, f32))> for LOR {
    fn from(((x1,y1,z1), (x2,y2,z2)): ((f32, f32, f32), (f32, f32, f32))) -> Self {
        Self { p1: Point::new(x1,y1,z1), p2: Point::new(x2,y2,z2), dt: 0.0  }
    }
}
