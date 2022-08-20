pub use crate::system_matrix::{LOR, find_tof_peak, find_entry_point, voxel_size, first_boundaries};

use ncollide3d as nc;
use units::todo::Lengthf32;

pub use geometry::{Vector, RatioPoint, RatioVec};

pub type Vectorf32 = nc::math::Vector<Lengthf32>;
pub type Pointf32  = nc::math::Point <Lengthf32>;
pub type Point     = geometry::Point;

pub use crate::index::{BoxDim_u, Index1_u, Index3_u, Index3Weightf32, LengthI, LengthU};

pub type BoundPair<T> = (std::ops::Bound<T>, std::ops::Bound<T>);
