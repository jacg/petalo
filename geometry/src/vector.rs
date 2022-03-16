use std::ops::{Index, Mul};
use uom::si::f32::{Length, Ratio};

//use crate::Point;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector {
    pub x: Length,
    pub y: Length,
    pub z: Length,
}



impl Mul<Ratio> for Vector {
    type Output = Self;
    fn mul(self, rhs: Ratio) -> Self::Output {
        Vector {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<Idx> Index<Idx> for Vector {
    type Output = Length;
    fn index(&self, _index: Idx) -> &Self::Output {
        todo!()
    }
}

impl Iterator for Vector {
    type Item = Length;

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

impl Vector
where
    Length: Mul,
{

    pub fn new(x: Length, y: Length, z: Length) -> Self { Self { x, y, z } }

    pub fn magnitude(&self) -> Length {
        // let Self { x, y, z } = self;
        // (x*x + y*y + z*z).sqrt()
        todo!()
    }

    pub fn argmin(self) -> (usize, Length) { todo!() }
    pub fn norm(self) -> Length { todo!() }
    pub fn normalize(self) -> Self { todo!() }

}
