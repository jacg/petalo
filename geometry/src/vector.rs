use std::ops::{Index, Mul};

//use crate::Point;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Vector<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}



impl<T> Mul<T> for Vector<T>
where
    T: Mul<Output = T> + Copy,
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        Vector {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T, Idx> Index<Idx> for Vector<T> {
    type Output = T;
    fn index(&self, _index: Idx) -> &Self::Output {
        todo!()
    }
}

impl<T> Iterator for Vector<T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

impl<T> Vector<T>
where
    T: Mul,
{

    pub fn new(x: T, y: T, z: T) -> Self { Self { x, y, z } }

    pub fn magnitude(&self) -> T {
        // let Self { x, y, z } = self;
        // (x*x + y*y + z*z).sqrt()
        todo!()
    }

    pub fn argmin(self) -> (usize, T) { todo!() }
    pub fn norm(self) -> T { todo!() }
    pub fn normalize(self) -> Self { todo!() }

}
