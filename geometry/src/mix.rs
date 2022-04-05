//! Utilities used in the transition to uom-aware types
//!
//! Everything in this module should probably be removed after the transition is
//! complete.

use crate::{Length, Point, Vector, RatioVec};


impl From<ncollide3d::math::Point<f32>> for Point {
    fn from(p: ncollide3d::math::Point<f32>) -> Self {
        use crate::uom::mm;
        let x: Length = mm(p.x);
        let y: Length = mm(p.y);
        let z: Length = mm(p.z);
        Self::new(x, y, z)
    }
}

impl From<ncollide3d::math::Vector<f32>> for Vector {
    fn from(v: ncollide3d::math::Vector<f32>) -> Self {
        use crate::uom::mm;
        let x: Length = mm(v.x);
        let y: Length = mm(v.y);
        let z: Length = mm(v.z);
        Self::new(x, y, z)
    }
}

impl From<Point> for ncollide3d::math::Point<f32> {
    fn from(p: Point) -> Self {
        use uom::si::length::millimeter as mm;
        let x = p.x.get::<mm>();
        let y = p.y.get::<mm>();
        let z = p.z.get::<mm>();
        Self::new(x, y, z)
    }
}

impl From<Vector> for ncollide3d::math::Vector<f32> {
    fn from(v: Vector) -> Self {
        use uom::si::length::millimeter as mm;
        let x = v.x.get::<mm>();
        let y = v.y.get::<mm>();
        let z = v.z.get::<mm>();
        Self::new(x, y, z)
    }
}


impl From<&ncollide3d::math::Point<f32>> for Point {
    fn from(p: &ncollide3d::math::Point<f32>) -> Self {
        use crate::uom::mm;
        let x: Length = mm(p.x);
        let y: Length = mm(p.y);
        let z: Length = mm(p.z);
        Self::new(x, y, z)
    }
}

impl From<&ncollide3d::math::Vector<f32>> for Vector {
    fn from(v: &ncollide3d::math::Vector<f32>) -> Self {
        use crate::uom::mm;
        let x: Length = mm(v.x);
        let y: Length = mm(v.y);
        let z: Length = mm(v.z);
        Self::new(x, y, z)
    }
}

impl From<&Point> for ncollide3d::math::Point<f32> {
    fn from(p: &Point) -> Self {
        use uom::si::length::millimeter as mm;
        let x = p.x.get::<mm>();
        let y = p.y.get::<mm>();
        let z = p.z.get::<mm>();
        Self::new(x, y, z)
    }
}

impl From<&Vector> for ncollide3d::math::Vector<f32> {
    fn from(v: &Vector) -> Self {
        use uom::si::length::millimeter as mm;
        let x = v.x.get::<mm>();
        let y = v.y.get::<mm>();
        let z = v.z.get::<mm>();
        Self::new(x, y, z)
    }
}


impl From<RatioVec> for ncollide3d::math::Vector<f32> {
    fn from(v: RatioVec) -> Self {
        use uom::si::ratio::ratio;
        let x = v.x.get::<ratio>();
        let y = v.y.get::<ratio>();
        let z = v.z.get::<ratio>();
        Self::new(x, y, z)
    }
}

impl From<&RatioVec> for ncollide3d::math::Vector<f32> {
    fn from(v: &RatioVec) -> Self {
        use uom::si::ratio::ratio;
        let x = v.x.get::<ratio>();
        let y = v.y.get::<ratio>();
        let z = v.z.get::<ratio>();
        Self::new(x, y, z)
    }
}

impl From<&ncollide3d::math::Vector<f32>> for RatioVec {
    fn from(v: &ncollide3d::math::Vector<f32>) -> Self {
        use crate::uom::ratio;
        let x = ratio(v.x);
        let y = ratio(v.y);
        let z = ratio(v.z);
        Self{x, y, z}
    }
}
