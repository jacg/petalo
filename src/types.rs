pub type Length = f32;
pub type Time   = Length;
pub type Weight = Length;
pub type Ratio  = Length;

pub type Index1 = usize;
pub type Index3 = [usize; 3];
pub type BoxDim = [usize; 3];

pub type Index3Weight = (Index3, Weight);

use ncollide3d as nc;
pub type Vector = nc::math ::Vector<Length>;
pub type Point  = nc::math ::Point <Length>;

// TODO: doesn't really belong in `types` ...
#[allow(clippy::excessive_precision)] // Stick to official definition of c
pub const C: Length = 0.299_792_458; // mm / ps

pub const TWOPI: Length = std::f32::consts::TAU as Length;
