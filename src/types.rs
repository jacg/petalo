use ncollide3d as nc;

pub use geometry::uom::uomcrate as guomc;
pub use guomc::si::{ISQ, SI, Quantity};
pub use guomc::typenum::{Z0, N1};
use geometry::in_base_unit;

pub type Lengthf32  = f32;
pub type UomLength  = geometry::Length;
pub type UomLengthI = geometry::uom::uomcrate::si::i32  ::Length;
pub type UomLengthU = geometry::uom::uomcrate::si::usize::Length;

pub type UomPerLength = geometry::PerLength;

pub type Timef32 = f32;
pub type UomTime = geometry::Time;

pub type UomVelocity = geometry::Velocity;

pub type Weightf32 = f32;  // TODO uom Weight

pub type Ratiof32 = f32;
pub type UomRatio = geometry::Ratio;

pub type Anglef32  = f32; // TODO uom Angle
pub type Energyf32 = f32; // TODO uom Energy
pub type Chargef32 = f32; // TODO uom Charge

pub type Intensityf32 = f32; // TODO uom Intensity

pub type Vectorf32 = nc::math::Vector<Lengthf32>;
pub type UomVector = geometry::Vector;

pub type    Point = nc::math::Point <Lengthf32>;
pub type UomPoint = geometry::Point;

pub type Index1 = usize;
pub type Index3 = [usize; 3];
pub type    BoxDim = [usize; 3];
pub type UomBoxDim = [UomLengthU; 3];

pub type Index3Weight = (Index3, Weightf32);

pub type BoundPair<T> = (std::ops::Bound<T>, std::ops::Bound<T>);

// TODO: doesn't really belong in `types` ...
#[allow(clippy::excessive_precision)] // Stick to official definition of c
pub const     C: Lengthf32 =               0.299_792_458; // mm / ps
pub const UOM_C: UomVelocity = in_base_unit!(299_792_458.0);

#[inline] pub fn ps_to_mm(dt: Timef32) -> Lengthf32 { dt * C }
#[inline] pub fn mm_to_ps(dx: Lengthf32) -> Timef32 { dx / C }

#[inline] pub fn ns_to_mm(dt: Timef32) -> Lengthf32 { ps_to_mm(dt) * 1000.0 }
#[inline] pub fn mm_to_ns(dx: Lengthf32) -> Timef32 { mm_to_ps(dx) / 1000.0  }

#[inline] pub fn ns_to_ps(dt: Timef32) -> Timef32 { dt * 1000.0 }


#[cfg(test)]
mod test_conversions {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    // Random test values
    const T: Timef32 = 3.5;
    const X: Lengthf32 = 1.2;

    #[test] fn human_a() { assert_approx_eq!(ns_to_mm(1.0), 300.0, 0.3  ); }
    #[test] fn human_b() { assert_approx_eq!(ps_to_mm(1.0),   0.3, 0.001); }

    #[test] fn relative_size_a() { assert_approx_eq!(ns_to_mm(T), 1000.0 * ps_to_mm(T)); }
    #[test] fn relative_size_b() { assert_approx_eq!(mm_to_ps(X), 1000.0 * mm_to_ns(X)); }

    #[test] fn roundtrip_a() { assert_approx_eq!(mm_to_ns(ns_to_mm(T)), T); }
    #[test] fn roundtrip_b() { assert_approx_eq!(mm_to_ps(ps_to_mm(T)), T); }

}

pub const     TWOPI: Lengthf32 = std::f32::consts::TAU as Lengthf32;
pub const UOM_TWOPI: UomRatio = in_base_unit!(std::f32::consts::TAU);
