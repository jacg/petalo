mod point;
mod vector;
pub mod units;

pub use point::{Point, RatioPoint};
pub use vector::{Vector, RatioVec};
pub use units::Quantity;
pub use units::{Angle, TWOPI, Length, Time, Velocity, Ratio, PerLength, AreaPerMass};

pub mod mix;

// Make the version of `uom` that is used here, accessible to other crates in
// the workspace. The problem is that the versions of `uom` declared as
// dependencies in differenc crates in the workspace, can diverge, and then we
// get annoying compilation errors. So, for now, we agree to use *this* version
// everywhere in the production code. The other version is used only in the
// `uom` example, which I want to leave in it's current prominent place, rather
// than moving it into this crate.
// TODO: this shouldn't be necessary any more once
// https://github.com/rust-lang/cargo/issues/8415 is stabilized.
pub use uom;
