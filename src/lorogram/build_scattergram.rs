use crate::Length;
use crate::lorogram::{Scattergram, axis_r, axis_phi, axis_z, axis_dz};
use ndhistogram::ndhistogram;
use geometry::uom::mm;


pub struct BuildScattergram {
    phi_bins: Option<usize>,
    r_bins:   Option<usize>,    z_bins: Option<usize>,
    dz_bins:  Option<usize>,    r_max:  Option<Length>,
    z_length: Option<Length>,  dz_max:  Option<Length>,
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
            dz_bins: None, dz_max: None,
        }
    }

    pub fn phi_bins(mut self, n: usize) -> Self { self.phi_bins = Some(n); self }

    pub fn r_bins(mut self, n: usize) -> Self {
        self.r_bins = Some(n);
        self.r_max.get_or_insert(mm(30.0));
        self
    }

    pub fn z_bins(mut self, n: usize) -> Self {
        self.z_bins = Some(n);
        self.z_length.get_or_insert(mm(200.0));
        self
    }

    pub fn dz_bins(mut self, n: usize) -> Self {
        self.dz_bins = Some(n);
        self.dz_max.get_or_insert(mm(1000.0));
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

    pub fn dz_max(mut self, z: Length) -> Self {
        self.dz_max = Some(z);
        self.dz_bins.get_or_insert(50);
        self
    }

    pub fn build(self) -> Option<Scattergram> {
        let phi = self.phi_bins;
        let r   = self.  r_bins.map(|n_bins| (n_bins, self. r_max  .unwrap()));
        let z   = self.  z_bins.map(|n_bins| (n_bins, self.z_length.unwrap()));
        let dz  = self. dz_bins.map(|n_bins| (n_bins, self.dz_max  .unwrap()));
        use Option::Some as S;
        let zax = |z_bins, z_length: Length| axis_z(z_bins, -z_length / 2.0, z_length / 2.0);
        match (r, z, phi, dz) {
            (None       , None       , None, None      ) => None,
            (None       , None       , None, S((db,dm))) => axes!(                                        axis_dz(db,dm)),
            (None       , None       , S(p), None      ) => axes!(                           axis_phi(p)                ),
            (None       , None       , S(p), S((db,dm))) => axes!(                           axis_phi(p), axis_dz(db,dm)),
            (None       , S((zb, zl)), None, None      ) => axes!(               zax(zb,zl)                             ),
            (None       , S((zb, zl)), None, S((db,dm))) => axes!(               zax(zb,zl)             , axis_dz(db,dm)),
            (None       , S((zb, zl)), S(p), None      ) => axes!(               zax(zb,zl), axis_phi(p)                ),
            (None       , S((zb, zl)), S(p), S((db,dm))) => axes!(               zax(zb,zl), axis_phi(p), axis_dz(db,dm)),
            (S((rb, rm)), None       , None, None      ) => axes!(axis_r(rb,rm)                                         ),
            (S((rb, rm)), None       , None, S((db,dm))) => axes!(axis_r(rb,rm)                         , axis_dz(db,dm)),
            (S((rb, rm)), None       , S(p), None      ) => axes!(axis_r(rb,rm),             axis_phi(p)                ),
            (S((rb, rm)), None       , S(p), S((db,dm))) => axes!(axis_r(rb,rm),             axis_phi(p), axis_dz(db,dm)),
            (S((rb, rm)), S((zb, zl)), None, None      ) => axes!(axis_r(rb,rm), zax(zb,zl)                             ),
            (S((rb, rm)), S((zb, zl)), None, S((db,dm))) => axes!(axis_r(rb,rm), zax(zb,zl)             , axis_dz(db,dm)),
            (S((rb, rm)), S((zb, zl)), S(p), None      ) => axes!(axis_r(rb,rm), zax(zb,zl), axis_phi(p)                ),
            (S((rb, rm)), S((zb, zl)), S(p), S((db,dm))) => axes!(axis_r(rb,rm), zax(zb,zl), axis_phi(p), axis_dz(db,dm)),
        }

    }

}
