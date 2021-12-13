pub type Length = f32;
pub type Time   = Length;
pub type Weight = Length;
pub type Ratio  = Length;
pub type Energy = Length;
pub type Charge = Length;
pub type Intensity = Length;

pub type Index1 = usize;
pub type Index3 = [usize; 3];
pub type BoxDim = [usize; 3];

pub type Index3Weight = (Index3, Weight);

use ncollide3d as nc;
pub type Vector = nc::math ::Vector<Length>;
pub type Point  = nc::math ::Point <Length>;

pub type BoundPair<T> = (std::ops::Bound<T>, std::ops::Bound<T>);

// TODO: doesn't really belong in `types` ...
#[allow(clippy::excessive_precision)] // Stick to official definition of c
pub const C: Length = 0.299_792_458; // mm / ps

#[inline] pub fn ps_to_mm(dt: Time) -> Length { dt * C }
#[inline] pub fn mm_to_ps(dx: Length) -> Time { dx / C }

#[inline] pub fn ns_to_mm(dt: Time) -> Length { ps_to_mm(dt) * 1000.0 }
#[inline] pub fn mm_to_ns(dx: Length) -> Time { mm_to_ps(dx) / 1000.0  }

#[inline] pub fn ns_to_ps(dt: Time) -> Time { dt * 1000.0 }


#[cfg(test)]
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

pub const TWOPI: Length = std::f32::consts::TAU as Length;
