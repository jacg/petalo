pub use crate::system_matrix::{LOR, find_tof_peak, find_entry_point, voxel_size, first_boundaries};

use ncollide3d as nc;

pub use geometry::uom::uomcrate as guomc;
pub use guomc::si::Quantity;
pub use guomc::typenum::{Z0, N1};
use geometry::in_base_unit;

pub type Lengthf32  = f32;
pub use geometry::Length;

pub type PerLength = geometry::PerLength;

pub type Timef32 = f32;
pub use geometry::Time;

pub type Velocity = geometry::Velocity;

pub type Weightf32 = f32;  // TODO uom Weight

pub type Ratiof32 = f32;
pub use geometry::Ratio;

pub use geometry::{Angle, TWOPI};
pub type Anglef32  = f32; // TODO eliminate

pub type Energyf32 = f32; // TODO uom Energy
pub type Chargef32 = f32; // TODO uom Charge

pub type Intensityf32 = f32; // TODO uom Intensity

pub use geometry::RatioVec;
pub type Vectorf32 = nc::math::Vector<Lengthf32>;
pub use geometry::Vector;

pub use geometry::RatioPoint;
pub type Pointf32 = nc::math::Point <Lengthf32>;
pub type Point    = geometry::Point;

pub use crate::index::{BoxDim_u, Index1_u, Index3_u, Index3Weightf32, LengthI, LengthU};

pub type BoundPair<T> = (std::ops::Bound<T>, std::ops::Bound<T>);

pub const C: Velocity = in_base_unit!(299_792_458.0);

pub use geometry::AreaPerMass;
