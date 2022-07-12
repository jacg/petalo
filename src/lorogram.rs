mod build_scattergram;
pub use build_scattergram::*;

use ndhistogram::{axis::{Axis, Uniform, UniformCyclic as Cyclic}, Histogram, ndhistogram};

use crate::system_matrix::LOR;
use std::f32::consts::TAU;

use crate::Lengthf32;
use crate::{Angle, Length, Point, Time, Ratio};
use geometry::units::{mm, mm_, ps_, ratio, radian_, turn};
use geometry::uom::ConstZero;


/// Distinguish between true, scatter and random prompt signals
pub enum Prompt { True, Scatter, Random }

pub struct Scattergram {
    trues   : Lorogram,
    scatters: Lorogram,
}

// TODO: consider adding a frozen version of the scattergram. This one needs two
// (trues, scatters) (or three (randoms)) histograms in order to accumulate
// data, but once all data have been collected, we only want the ratios of the
// bin values, so keeping multiple separate histograms is a waste of memory, and
// computing the ratios repeatedly on the fly is a waste of time.
impl Scattergram {

    pub fn new(
        bins_phi: usize,
        bins_z  : usize, len_z : Length,
        bins_dz : usize, len_dz: Length,
        bins_r  : usize, max_r : Length,
        bins_dt : usize, max_dt: Time
    ) -> Self {
        let max_z = len_z / 2.0;
        let trues = Lorogram(ndhistogram!(
            LorAxisPhi::new(bins_phi),
            LorAxisZ::new(bins_z, -max_z, max_z),
            LorAxisDz::new(bins_dz, len_dz),
            LorAxisR::new(bins_r, max_r),
            LorAxisT::new(bins_dt, max_dt);
            usize
        ));
        // TODO: Can we clone `trues`?
        let scatters = Lorogram(ndhistogram!(
            LorAxisPhi::new(bins_phi),
            LorAxisZ::new(bins_z, -max_z, max_z),
            LorAxisDz::new(bins_dz, len_dz),
            LorAxisR::new(bins_r, max_r),
            LorAxisT::new(bins_dt, max_dt);
            usize
        ));
        Self { trues, scatters }
    }

    pub fn fill(&mut self, kind: Prompt, lor: &LOR) {
        match kind {
            Prompt::True    => self.trues.   fill(lor),
            Prompt::Scatter => self.scatters.fill(lor),
            Prompt::Random  => panic!("Not expecting any random events yet."),
        }
    }

    /// Multiplicative contribution of scatters to trues, in nearby LORs.
    ///
    /// `(scatters + trues) / trues`
    pub fn value(&self, lor: &LOR) -> Ratio {
        let trues = self.trues.value(lor);
        ratio(if trues > 0 {
            let scatters: f32 = self.scatters.value(lor) as f32;
            let trues = trues as f32;
            (scatters + trues) / trues
        } else { f32::MAX })
    }

    pub fn triplet(&self, lor: &LOR) -> (Ratio, f32, f32) {
        let trues = self.trues.value(lor);
        if trues > 0 {
            let scatters: f32 = self.scatters.value(lor) as f32;
            let trues = trues as f32;
            (ratio((scatters + trues) / trues), trues, scatters)
        } else { (ratio(1.0), 0.0, self.scatters.value(lor) as f32) }
    }
}
// --------------------------------------------------------------------------------
pub struct MappedAxis<T,A>
where
    A: Axis,
{
    axis: A,
    map: Box<dyn Fn(&T) -> A::Coordinate + Sync>,
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
pub struct LorAxisPhi {
    axis: Cyclic<Lengthf32>,
}
impl LorAxisPhi {
    fn new(nbins: usize) -> Self {
        Self { axis: Cyclic::new(nbins, 0.0, TAU), }
    }
}
impl Axis for LorAxisPhi {
    type Coordinate = LOR;
    type BinInterval = <Uniform<Lengthf32> as Axis>::BinInterval;
    fn index(&self, lor: &LOR) -> Option<usize> { self.axis.index(&radian_(phi(lor))) }
    fn num_bins(&self) -> usize { self.axis.num_bins() }
    fn bin(&self, index: usize) -> Option<Self::BinInterval> { self.axis.bin(index) }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pub struct LorAxisZ {
    axis: Uniform<Lengthf32>,
}
impl LorAxisZ {
    fn new(nbins: usize, min: Length, max: Length) -> Self {
        Self { axis: Uniform::new(nbins, mm_(min), mm_(max)), }
    }
}
impl Axis for LorAxisZ {
    type Coordinate = LOR;
    type BinInterval = <Uniform<Lengthf32> as Axis>::BinInterval;
    fn index(&self, lor: &LOR) -> Option<usize> { self.axis.index(&mm_(z_of_midpoint(lor))) }
    fn num_bins(&self) -> usize { self.axis.num_bins() }
    fn bin(&self, index: usize) -> Option<Self::BinInterval> { self.axis.bin(index) }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pub struct LorAxisDz {
    axis: Uniform<Lengthf32>,
}
impl LorAxisDz {
    fn new(nbins: usize, len: Length) -> Self {
        Self { axis: Uniform::new(nbins, 0.0, mm_(len)), }
    }
}
impl Axis for LorAxisDz {
    type Coordinate = LOR;
    type BinInterval = <Uniform<Lengthf32> as Axis>::BinInterval;
    fn index(&self, lor: &LOR) -> Option<usize> { self.axis.index(&mm_(delta_z(lor))) }
    fn num_bins(&self) -> usize { self.axis.num_bins() }
    fn bin(&self, index: usize) -> Option<Self::BinInterval> { self.axis.bin(index) }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pub struct LorAxisR {
    axis: Uniform<Lengthf32>,
}
impl LorAxisR {
    fn new(nbins: usize, len: Length) -> Self {
        Self { axis: Uniform::new(nbins, 0.0, mm_(len)), }
    }
}
impl Axis for LorAxisR {
    type Coordinate = LOR;
    type BinInterval = <Uniform<Lengthf32> as Axis>::BinInterval;
    fn index(&self, lor: &LOR) -> Option<usize> { self.axis.index(&mm_(distance_from_z_axis(lor))) }
    fn num_bins(&self) -> usize { self.axis.num_bins() }
    fn bin(&self, index: usize) -> Option<Self::BinInterval> { self.axis.bin(index) }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pub struct LorAxisT {
    axis: Uniform<Lengthf32>,
}
impl LorAxisT {
    fn new(nbins: usize, max: Time) -> Self {
        Self { axis: Uniform::new(nbins, ps_(-max), ps_(max)), }
    }
}
impl Axis for LorAxisT {
    type Coordinate = LOR;
    type BinInterval = <Uniform<Lengthf32> as Axis>::BinInterval;
    fn index(&self, lor: &LOR) -> Option<usize> { self.axis.index(&mm_(distance_from_z_axis(lor))) }
    fn num_bins(&self) -> usize { self.axis.num_bins() }
    fn bin(&self, index: usize) -> Option<Self::BinInterval> { self.axis.bin(index) }
}
// ================================================================================
struct Lorogram(ndhistogram::HistND<(LorAxisPhi, LorAxisZ, LorAxisDz, LorAxisR, LorAxisT), usize>);

impl Lorogram {
    pub fn fill (&mut self, lor: &LOR)          {  self.0.fill (&(*lor, *lor, *lor, *lor, *lor))               }
    pub fn value(&    self, lor: &LOR) -> usize { *self.0.value(&(*lor, *lor, *lor, *lor, *lor)).unwrap_or(&0) }
}


pub type LorAxU = MappedAxis<LOR, Uniform<Lengthf32>>;
pub type LorAxC = MappedAxis<LOR, Cyclic <Lengthf32>>;

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
    // TODO this repeats the work done in distance_from_z_axis. Can this be
    // optimized out, once we settle on a less flexible scattergram?
    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;
    let x1 = p1.x;
    let y1 = p1.y;
    let r = (dx * y1 - dy * x1) / (dx*dx + dy*dy).sqrt();
    let phi = phi_of_x_y(dx, dy);
    if r < mm(0.0) { phi + turn(0.5) }
    else           { phi             }

}

fn phi_of_x_y(x: Length, y: Length) -> Angle { y.atan2(x) }


pub fn mk_lor(((x1,y1,z1), (x2,y2,z2)): ((f32, f32, f32), (f32, f32, f32))) -> LOR {
    let (x1, y1, z1, x2, y2, z2) = (mm(x1), mm(y1), mm(z1), mm(x2), mm(y2), mm(z2));
    LOR { p1: Point::new(x1,y1,z1), p2: Point::new(x2,y2,z2), dt: Time::ZERO, additive_correction: ratio(1.0) }
}
