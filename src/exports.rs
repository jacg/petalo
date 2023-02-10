pub use crate::{
    index::{BoxDim_u, Index1_u, Index3_u, LengthI, LengthU},
    lor::LOR,
};

use ncollide3d as nc;
use units::todo::Lengthf32;

pub use geometry::{Vector, RatioPoint, RatioVec};

pub type Vectorf32 = nc::math::Vector<Lengthf32>;
pub type Pointf32  = nc::math::Point <Lengthf32>;
pub type Point     = geometry::Point;


pub type BoundPair<T> = (std::ops::Bound<T>, std::ops::Bound<T>);
