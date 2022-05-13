use crate::Length;
use crate::lorogram::{Scattergram, axis_r, axis_phi};
use ndhistogram::ndhistogram;
use geometry::uom::mm;


pub struct BuildScattergram {
    phi_bins: Option<usize>,
    r_bins: Option<usize>,
    r_max: Option<Length>,
}

impl BuildScattergram {

    pub fn new() -> Self { Self { phi_bins: None, r_bins: None, r_max: None } }

    pub fn phi_bins(mut self, n: usize)  -> Self { self.phi_bins = Some(n); self }

    pub fn r_bins(mut self, n: usize)  -> Self {
        self.r_bins = Some(n);
        self.r_max.get_or_insert(mm(30.0));
        self
    }
    pub fn r_max(mut self, r: Length) -> Self {
        self.r_max = Some(r);
        self.r_bins.get_or_insert(50);
        self
    }

    pub fn build(self) -> Scattergram {
        let r   = self.  r_bins.map(|n_bins| (n_bins, self.r_max.unwrap()));
        let phi = self.phi_bins;
        let f = Scattergram::new;
        if let Some((r_bins, r_max)) = r {
            if let Some(n) = phi { f(&|| Box::new(ndhistogram!(axis_r(r_bins, r_max), axis_phi(n); usize))) }
            else                 { f(&|| Box::new(ndhistogram!(axis_r(r_bins, r_max)             ; usize))) }
        } else {
            if let Some(n) = phi { f(&|| Box::new(ndhistogram!(                       axis_phi(n); usize))) }
            else                   { panic!("Histogram without axes") }
        }
    }

}
