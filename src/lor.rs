use units::{mm_, ns_, C, Length, Ratio, Time};
use crate::Point;


/// Line Of Response.
///
/// 2 spacetime vectors indicating the positions and times of coincident
/// detector element activations
#[derive(Clone, Copy, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct LOR {
    pub p1: Point,
    pub p2: Point,
    pub dt: Time,
    /// Scatter and random corrections, which appear as an additive contribution
    /// to the sinogram bin in the MLEM forward projection. In order to make it
    /// compatible with a single LOR (rather than many LORs in a sinogram bin)
    /// it is expressed here as a *multiplicative* factor.
    pub additive_correction: Ratio,
}

impl LOR {
    pub fn new(t1: Time, t2: Time, p1: Point, p2: Point, additive_correction: Ratio) -> Self {
        Self { p1, p2, dt: t2 - t1, additive_correction }
    }

    pub fn from_components((t1, t2): (Time, Time),
                           (x1, y1, z1): (Length, Length, Length),
                           (x2, y2, z2): (Length, Length, Length),
                           additive_correction: Ratio
                          ) -> Self
    {
        Self::new(t1, t2, Point::new(x1,y1,z1), Point::new(x2,y2,z2), additive_correction)
    }

}

use core::fmt;
impl fmt::Display for LOR {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let (p, q) = (self.p1, self.p2);
        write!(f, "<LOR ({:8.2} {:8.2} {:8.2}) ({:8.2} {:8.2} {:8.2}) {:7.2}ns {:7.2}mm /{:7.2} >",
               mm_(p.x), mm_(p.y), mm_(p.z),
               mm_(q.x), mm_(q.y), mm_(q.z),
               ns_(self.dt), mm_(self.dt * C) / 2.0,
               mm_((p-q).norm())
        )
    }
}
