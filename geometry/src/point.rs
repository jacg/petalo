use std::ops::{AddAssign, Index, Sub};
use uom::si::f32::Length;
use crate::Vector;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point {
    pub x: Length,
    pub y: Length,
    pub z: Length,
}

impl Point {
    pub fn new(x: Length, y: Length, z: Length) -> Self { Self { x, y, z } }
    pub fn component_div(&mut self, _other: &Self) -> Self { todo!() }
}

impl Sub for Point {
    type Output = Vector;
    fn sub(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Sub for &Point {
    type Output = Vector;
    fn sub(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl AddAssign<Vector> for Point {
    fn add_assign(&mut self, _rhs: Vector) {
        todo!()
    }
}

impl<Idx> Index<Idx> for Point {
    type Output = Length;
    fn index(&self, _index: Idx) -> &Self::Output {
        todo!()
    }
}
