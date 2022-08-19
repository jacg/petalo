use units::{Length, Time, mm, ps};
use crate::lorogram::Scattergram;

pub struct BuildScattergram {
    phi_bins: Option<usize>,
    r_bins  : Option<usize>, r_max   : Option<Length>,
    z_bins  : Option<usize>, z_length: Option<Length>,
    dz_bins : Option<usize>, dz_max  : Option<Length>,
    dt_bins : Option<usize>, dt_max  : Option<Time>,
//
// NOTE: Fine-grained bins seem to give bad reconstructed images: perhaps too
// low statistics. If this is the case, then interpolation in Scattergram::value
// could be important.
}

impl From<&crate::config::mlem::Scatter> for Option<Scattergram> {
    fn from(config: &crate::config::mlem::Scatter) -> Self {
        use crate::config::mlem::{Bins, BinsMax, BinsLength};

        let mut builder = BuildScattergram::new();

        if let Some(Bins { bins }) = config.phi {
            builder = builder.phi_bins(bins);
        }
        if let Some(BinsMax { bins, max }) = config.r {
            builder = builder.r_bins(bins);
            builder = builder.r_max (max );
        }
        if let Some(BinsMax { bins, max }) = config.dz {
            builder = builder.dz_bins(bins);
            builder = builder.dz_max (max );
        }
        if let Some(BinsMax { bins, max }) = config.dt {
            builder = builder.dt_bins(bins);
            builder = builder.dt_max (max );
        }
        if let Some(BinsLength { bins, length }) = config.z {
            builder = builder.z_bins  (bins);
            builder = builder.z_length(length);
        }
        builder.build()
    }
}

const DEFAULT_NUMBER_OF_BINS: usize = 30;

impl BuildScattergram {

    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self {
            phi_bins: None,
            r_bins  : None, r_max   : None,
            z_bins  : None, z_length: None,
            dz_bins : None, dz_max  : None,
            dt_bins : None, dt_max  : None,
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

    pub fn dt_bins(mut self, n: usize) -> Self {
        self.dt_bins = Some(n);
        self.dt_max.get_or_insert(ps(1700.0)); // Δ-TOF ~50cm off-centre
        self
    }

    pub fn r_max(mut self, r: Length) -> Self {
        self.r_max = Some(r);
        self.r_bins.get_or_insert(DEFAULT_NUMBER_OF_BINS);
        self
    }

    pub fn z_length(mut self, z: Length) -> Self {
        self.z_length = Some(z);
        self.z_bins.get_or_insert(DEFAULT_NUMBER_OF_BINS);
        self
    }

    pub fn dz_max(mut self, z: Length) -> Self {
        self.dz_max = Some(z);
        self.dz_bins.get_or_insert(DEFAULT_NUMBER_OF_BINS);
        self
    }

    pub fn dt_max(mut self, dt: Time) -> Self {
        self.dt_max = Some(dt);
        self.dt_bins.get_or_insert(DEFAULT_NUMBER_OF_BINS);
        self
    }

    pub fn build(self) -> Option<Scattergram> {
        let phi = self.phi_bins;
        let r   = self.  r_bins.map(|n_bins| (n_bins, self. r_max  .unwrap()));
        let z   = self.  z_bins.map(|n_bins| (n_bins, self.z_length.unwrap()));
        let dz  = self. dz_bins.map(|n_bins| (n_bins, self.dz_max  .unwrap()));
        let dt  = self. dt_bins.map(|n_bins| (n_bins, self.dt_max  .unwrap()));
        // TODO: tepmorary hack for default values
        let  l  = mm(66666666.6);
        let  t  = ps(66666666.6);
        use Option::Some as S;
        let new = Scattergram::new;
        match (dt, r, z, phi, dz) {
            (None       , None       , None       , None, None      ) => None,
            (None       , None       , None       , None, S((db,dm))) => Some(new(1,  1,  l, db, dm,  1,  l,  1,  t)),
            (None       , None       , None       , S(p), None      ) => Some(new(p,  1,  l,  1,  l,  1,  l,  1,  t)),
            (None       , None       , None       , S(p), S((db,dm))) => Some(new(p,  1,  l, db, dm,  1,  l,  1,  t)),
            (None       , None       , S((zb, zl)), None, None      ) => Some(new(1, zb, zl,  1,  l,  1,  l,  1,  t)),
            (None       , None       , S((zb, zl)), None, S((db,dm))) => Some(new(1, zb, zl, db, dm,  1,  l,  1,  t)),
            (None       , None       , S((zb, zl)), S(p), None      ) => Some(new(p, zb, zl,  1,  l,  1,  l,  1,  t)),
            (None       , None       , S((zb, zl)), S(p), S((db,dm))) => Some(new(p, zb, zl, db, dm,  1,  l,  1,  t)),
            (None       , S((rb, rm)), None       , None, None      ) => Some(new(1,  1,  l,  1,  l, rb, rm,  1,  t)),
            (None       , S((rb, rm)), None       , None, S((db,dm))) => Some(new(1,  1,  l, db, dm, rb, rm,  1,  t)),
            (None       , S((rb, rm)), None       , S(p), None      ) => Some(new(p,  1,  l,  1,  l, rb, rm,  1,  t)),
            (None       , S((rb, rm)), None       , S(p), S((db,dm))) => Some(new(p,  1,  l, db, dm, rb, rm,  1,  t)),
            (None       , S((rb, rm)), S((zb, zl)), None, None      ) => Some(new(1, zb, zl,  1,  l, rb, rm,  1,  t)),
            (None       , S((rb, rm)), S((zb, zl)), None, S((db,dm))) => Some(new(1, zb, zl, db, dm, rb, rm,  1,  t)),
            (None       , S((rb, rm)), S((zb, zl)), S(p), None      ) => Some(new(p, zb, zl,  1,  l, rb, rm,  1,  t)),
            (None       , S((rb, rm)), S((zb, zl)), S(p), S((db,dm))) => Some(new(p, zb, zl, db, dm, rb, rm,  1,  t)),
            (S((tb, tm)), None       , None       , None, None      ) => Some(new(1,  1,  l,  1,  l,  1,  l, tb, tm)),
            (S((tb, tm)), None       , None       , None, S((db,dm))) => Some(new(1,  1,  l, db, dm,  1,  l, tb, tm)),
            (S((tb, tm)), None       , None       , S(p), None      ) => Some(new(p,  1,  l,  1,  l,  1,  l, tb, tm)),
            (S((tb, tm)), None       , None       , S(p), S((db,dm))) => Some(new(p,  1,  l, db, dm,  1,  l, tb, tm)),
            (S((tb, tm)), None       , S((zb, zl)), None, None      ) => Some(new(1, zb, zl,  1,  l,  1,  l, tb, tm)),
            (S((tb, tm)), None       , S((zb, zl)), None, S((db,dm))) => Some(new(1, zb, zl, db, dm,  1,  l, tb, tm)),
            (S((tb, tm)), None       , S((zb, zl)), S(p), None      ) => Some(new(p, zb, zl,  1,  l,  1,  l, tb, tm)),
            (S((tb, tm)), None       , S((zb, zl)), S(p), S((db,dm))) => Some(new(p, zb, zl, db, dm,  1,  l, tb, tm)),
            (S((tb, tm)), S((rb, rm)), None       , None, None      ) => Some(new(1,  1,  l,  1,  l, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), None       , None, S((db,dm))) => Some(new(1,  1,  l, db, dm, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), None       , S(p), None      ) => Some(new(p,  1,  l,  1,  l, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), None       , S(p), S((db,dm))) => Some(new(p,  1,  l, db, dm, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), S((zb, zl)), None, None      ) => Some(new(1, zb, zl,  1,  l, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), S((zb, zl)), None, S((db,dm))) => Some(new(1, zb, zl, db, dm, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), S((zb, zl)), S(p), None      ) => Some(new(p, zb, zl,  1,  l, rb, rm, tb, tm)),
            (S((tb, tm)), S((rb, rm)), S((zb, zl)), S(p), S((db,dm))) => Some(new(p, zb, zl, db, dm, rb, rm, tb, tm)),
        }

    }

}
