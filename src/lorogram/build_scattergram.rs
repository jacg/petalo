use crate::Length;
use crate::lorogram::{Scattergram, axis_r, axis_phi, axis_z};
use ndhistogram::ndhistogram;
use geometry::uom::mm;


pub struct BuildScattergram {
    phi_bins: Option<usize>,
    r_bins:   Option<usize>,
    r_max:    Option<Length>,
    z_bins:   Option<usize>,
    z_length: Option<Length>,
//
// NOTE: Fine-grained bins seem to give bad reconstructed images: perhaps too
// low statistics. If this is the case, then interpolation in Scattergram::value
// could be important.


}

macro_rules! axes {
    ($($axes:expr),+) => {
        Some(Scattergram::new(&|| Box::new(ndhistogram!($($axes),+; usize))))
    };
}

impl BuildScattergram {

    pub fn new() -> Self {
        Self {
            phi_bins: None,
            r_bins: None, r_max: None,
            z_bins: None, z_length: None,
        }
    }

    pub fn phi_bins(mut self, n: usize)  -> Self { self.phi_bins = Some(n); self }

    pub fn r_bins(mut self, n: usize)  -> Self {
        self.r_bins = Some(n);
        self.r_max.get_or_insert(mm(30.0));
        self
    }

    pub fn z_bins(mut self, n: usize)  -> Self {
        self.z_bins = Some(n);
        self.z_length.get_or_insert(mm(200.0));
        self
    }

    pub fn r_max(mut self, r: Length) -> Self {
        self.r_max = Some(r);
        self.r_bins.get_or_insert(50);
        self
    }

    pub fn z_length(mut self, z: Length) -> Self {
        self.z_length = Some(z);
        self.z_bins.get_or_insert(50);
        self
    }

    pub fn build(self) -> Option<Scattergram> {
        let r   = self.  r_bins.map(|n_bins| (n_bins, self.r_max   .unwrap()));
        let z   = self.  z_bins.map(|n_bins| (n_bins, self.z_length.unwrap()));
        let phi = self.phi_bins;
        use Option::Some as S;
        let zax = |z_bins, z_length: Length| axis_z(z_bins, -z_length / 2.0, z_length / 2.0);
        match (r, z, phi) {
            (None       , None       , None) => None,
            (None       , None       , S(p)) => axes!(                           axis_phi(p)),
            (None       , S((zb, zl)), None) => axes!(               zax(zb,zl)             ),
            (None       , S((zb, zl)), S(p)) => axes!(               zax(zb,zl), axis_phi(p)),
            (S((rb, rm)), None       , None) => axes!(axis_r(rb,rm)                         ),
            (S((rb, rm)), None       , S(p)) => axes!(axis_r(rb,rm),             axis_phi(p)),
            (S((rb, rm)), S((zb, zl)), None) => axes!(axis_r(rb,rm), zax(zb,zl)             ),
            (S((rb, rm)), S((zb, zl)), S(p)) => axes!(axis_r(rb,rm), zax(zb,zl), axis_phi(p)),
        }

    }

}
