//! Utilities for converting between uom-aware and uom-unaware types.
//!
//! This module should be removed once the transition to uom-awareness is
//! complete.

use crate::types::{Point, Vector};
use crate::types::{UomLength, UomPoint, UomVector};

pub fn uom_vector_to_ncollide(v: &UomVector) -> Vector {
    use geometry::uom::uomcrate::si::length::millimeter as mm;
    let x = v.x.get::<mm>();
    let y = v.y.get::<mm>();
    let z = v.z.get::<mm>();
    Vector::new(x, y, z)
}

pub fn uom_point_as_ncollide(p: &UomPoint) -> Point {
    use geometry::uom::uomcrate::si::length::millimeter as mm;
    let x = p.x.get::<mm>();
    let y = p.y.get::<mm>();
    let z = p.z.get::<mm>();
    Point::new(x, y, z)
}

pub fn uom_vector_from_ncollide(p: &Vector) -> UomVector {
    use geometry::uom::mm;
    let x: UomLength = mm(p.x);
    let y: UomLength = mm(p.y);
    let z: UomLength = mm(p.z);
    UomVector::new(x, y, z)
}

pub fn uom_point_from_ncollide(p: &Point) -> UomPoint {
    use geometry::uom::mm;
    let x: UomLength = mm(p.x);
    let y: UomLength = mm(p.y);
    let z: UomLength = mm(p.z);
    UomPoint::new(x, y, z)
}
