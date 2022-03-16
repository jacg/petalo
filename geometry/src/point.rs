use std::ops::{AddAssign, Index, Sub};
use crate::Vector;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Point<T> {
    pub fn new(x: T, y: T, z: T) -> Self { Self { x, y, z } }
    pub fn component_div(&mut self, _other: &Self) -> Self { todo!() }
}

impl<T> Sub for Point<T>
where
    T: Sub<Output = T>,
{
    type Output = Vector<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Sub for &Point<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Vector<T>;
    fn sub(self, rhs: Self) -> Self::Output {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> AddAssign<Vector<T>> for Point<T> {
    fn add_assign(&mut self, _rhs: Vector<T>) {
        todo!()
    }
}

impl<T, Idx> Index<Idx> for Point<T> {
    type Output = T;
    fn index(&self, _index: Idx) -> &Self::Output {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_minus_point() {
        let pl = Point::<i32>::new(10, 20, 30);
        let pr = Point::<i32>::new( 1,  2,  3);
        assert_eq!(pl - pr, Vector::<i32>::new(9, 18, 27));
    }
}
