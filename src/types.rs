#[cfg    (feature = "units") ] pub use geometry::uom::uomcrate as guomc;
#[cfg    (feature = "units") ] pub use guomc::si::{ISQ, SI, Quantity};
#[cfg    (feature = "units") ] pub use guomc::typenum::{Z0, N1};


#[cfg(not(feature = "units"))] pub type  Length = f32;
#[cfg    (feature = "units") ] pub type  Length = geometry::Length;
#[cfg    (feature = "units") ] pub type ILength = geometry::uom::uomcrate::si::i32  ::Length;
#[cfg    (feature = "units") ] pub type ULength = geometry::uom::uomcrate::si::usize::Length;

#[cfg(not(feature = "units"))] pub type PerLength = f32;
#[cfg    (feature = "units") ] pub type PerLength = Quantity<ISQ<N1, Z0, Z0, Z0, Z0, Z0, Z0>, SI<f32>, f32>;

#[cfg(not(feature = "units"))] pub type Time = Length;
#[cfg    (feature = "units") ] pub type Time = geometry::Time;

#[cfg    (feature = "units") ] pub type Velocity = geometry::Velocity;

#[cfg(not(feature = "units"))] pub type Weight = f32;
#[cfg    (feature = "units") ] pub type Weight = f32;  // TODO uom Weight

#[cfg(not(feature = "units"))] pub type Ratio = f32;
#[cfg    (feature = "units") ] pub type Ratio = geometry::Ratio;

#[cfg(not(feature = "units"))] pub type Angle = f32;
#[cfg    (feature = "units") ] pub type Angle = f32; // TODO uom Angl

#[cfg(not(feature = "units"))] pub type Energy = f32;
#[cfg    (feature = "units") ] pub type Energy = f32; //TODO uom Energy

#[cfg(not(feature = "units"))] pub type Charge = f32;
#[cfg    (feature = "units") ] pub type Charge = f32; // TODO uom Charge

#[cfg(not(feature = "units"))] pub type Intensity = f32;
#[cfg    (feature = "units") ] pub type Intensity = f32; // TODO uom Intensity

#[cfg(not(feature = "units"))] pub type Vector = nc::math::Vector<Length>;
#[cfg    (feature = "units") ] pub type Vector = geometry::Vector;

#[cfg(not(feature = "units"))] pub type Point  = nc::math::Point <Length>;
#[cfg    (feature = "units") ] pub type Point  = geometry::Point;

pub type Index1 = usize;
pub type Index3 = [usize; 3];
pub type BoxDim = [usize; 3];

pub type Index3Weight = (Index3, Weight);

#[cfg(not(feature = "units"))] use ncollide3d as nc;

pub type BoundPair<T> = (std::ops::Bound<T>, std::ops::Bound<T>);

// TODO: doesn't really belong in `types` ...
#[allow(clippy::excessive_precision)] // Stick to official definition of c
#[cfg(not(feature = "units"))] pub const C: Length = 0.299_792_458; // mm / ps
#[cfg    (feature = "units") ]
pub const C: Velocity = geometry::Quantity {
  dimension: std::marker::PhantomData,
  units: std::marker::PhantomData,
  value: 299_792_458.0, // How can I be sure that this is in m/s ?
};

#[inline] pub fn ps_to_mm(dt: Time) -> Length { dt * C }
#[inline] pub fn mm_to_ps(dx: Length) -> Time { dx / C }

#[inline] pub fn ns_to_mm(dt: Time) -> Length { ps_to_mm(dt) * 1000.0 }
#[inline] pub fn mm_to_ns(dx: Length) -> Time { mm_to_ps(dx) / 1000.0  }

#[inline] pub fn ns_to_ps(dt: Time) -> Time { dt * 1000.0 }


#[cfg(test)]
#[cfg(not(feature = "units"))]
mod test_conversions {
    use super::*;
    use assert_approx_eq::assert_approx_eq;

    // Random test values
    const T: Time = 3.5;
    const X: Length = 1.2;

    #[test] fn human_a() { assert_approx_eq!(ns_to_mm(1.0), 300.0, 0.3  ); }
    #[test] fn human_b() { assert_approx_eq!(ps_to_mm(1.0),   0.3, 0.001); }

    #[test] fn relative_size_a() { assert_approx_eq!(ns_to_mm(T), 1000.0 * ps_to_mm(T)); }
    #[test] fn relative_size_b() { assert_approx_eq!(mm_to_ps(X), 1000.0 * mm_to_ns(X)); }

    #[test] fn roundtrip_a() { assert_approx_eq!(mm_to_ns(ns_to_mm(T)), T); }
    #[test] fn roundtrip_b() { assert_approx_eq!(mm_to_ps(ps_to_mm(T)), T); }

}

#[cfg(not(feature = "units"))] pub const TWOPI: Length = std::f32::consts::TAU as Length;
#[cfg    (feature = "units") ]
pub const TWOPI: Ratio = geometry::Quantity {
  dimension: std::marker::PhantomData,
  units: std::marker::PhantomData,
  value: std::f32::consts::TAU,
};
