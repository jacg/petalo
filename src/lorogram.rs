mod build_scattergram;
pub use build_scattergram::*;

use ndhistogram::{axis::{Axis, Uniform, UniformCyclic as Cyclic},
                  Histogram, ndhistogram};

use crate::system_matrix::LOR;
use std::f32::consts::TAU;

use crate::Lengthf32;
use units::{Angle, Length, Time, Ratio};
use units::{mm, mm_, ps_, ratio, radian_, turn};


/// Distinguish between true, scatter and random prompt signals
pub enum Prompt { True, Scatter, Random }

#[derive(Clone)]
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

    #[allow(clippy::too_many_arguments)]
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
            LorAxisZ  ::new(bins_z, -max_z, max_z),
            LorAxisDz ::new(bins_dz, len_dz),
            LorAxisR  ::new(bins_r,  max_r ),
            LorAxisT  ::new(bins_dt, max_dt);
            usize
        ));
        // TODO: Can we clone `trues`?
        let scatters = trues.clone();
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


impl std::ops::AddAssign<&Scattergram> for Scattergram {
    fn add_assign(&mut self, rhs: &Self) {
        self.trues    += &rhs.trues;
        self.scatters += &rhs.scatters;
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
macro_rules! lor_axis {
    ($name:ident<$axis_kind:ident> {
        fn new($($parameter:tt: $type:ty),*) { axis::new( $($new_arg:expr),* )}
        index: $coord:ident -> $index:expr
    }) => {
        #[derive(Clone, PartialEq)]
        pub struct $name {
            axis: $axis_kind<Lengthf32>,
        }
        impl $name {
            fn new($($parameter: $type),*) -> Self {
                Self { axis: $axis_kind::new($($new_arg),*) }
            }
        }
        impl Axis for $name {
            type Coordinate = LOR;
            type BinInterval = <Uniform<Lengthf32> as Axis>::BinInterval;
            fn index(&self, $coord: &LOR) -> Option<usize> { self.axis.index(&$index) }
            fn num_bins(&self) -> usize { self.axis.num_bins() }
            fn bin(&self, index: usize) -> Option<Self::BinInterval> { self.axis.bin(index) }
        }
    };
}

lor_axis!{
    LorAxisPhi<Cyclic> {
        fn new(nbins: usize) { axis::new(nbins, 0.0, TAU) }
        index: lor -> radian_(phi(lor))
    }
}

lor_axis!{
    LorAxisZ<Uniform> {
        fn new(nbins: usize, min: Length, max: Length) { axis::new(nbins, mm_(min), mm_(max)) }
        index: lor -> mm_(z_of_midpoint(lor))
    }
}

lor_axis!{
    LorAxisDz<Uniform> {
        fn new(nbins: usize, len: Length) { axis::new(nbins, 0.0, mm_(len)) }
        index: lor -> mm_(delta_z(lor))
    }
}

lor_axis!{
    LorAxisR<Uniform> {
        fn new(nbins: usize, len: Length) { axis::new(nbins, 0.0, mm_(len)) }
        index: lor -> mm_(distance_from_z_axis(lor))
    }
}

lor_axis!{
    LorAxisT<Uniform> {
        fn new(nbins: usize, max: Time) { axis::new(nbins, ps_(-max), ps_(max)) }
        index: lor -> mm_(distance_from_z_axis(lor))
    }
}
// ================================================================================
#[derive(Clone)]
struct Lorogram(ndhistogram::HistND<(LorAxisPhi, LorAxisZ, LorAxisDz, LorAxisR, LorAxisT), usize>);

impl Lorogram {
    pub fn fill (&mut self, lor: &LOR)          {  self.0.fill (&(*lor, *lor, *lor, *lor, *lor))               }
    pub fn value(&    self, lor: &LOR) -> usize { *self.0.value(&(*lor, *lor, *lor, *lor, *lor)).unwrap_or(&0) }
}

impl std::ops::AddAssign<&Lorogram> for Lorogram {
    fn add_assign(&mut self, rhs: &Lorogram) {
        self.0 += &rhs.0;
    }
}


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
