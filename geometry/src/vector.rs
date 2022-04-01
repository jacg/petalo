use std::ops::{Index, Mul};
use crate::uom::mmps::f32::Length;

//use crate::Point;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vector {
    pub x: Length,
    pub y: Length,
    pub z: Length,
}



impl Mul<f32> for Vector {
    type Output = Self;
    fn mul(self, rhs: f32) -> Self::Output {
        Vector {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Index<usize> for Vector {
    type Output = Length;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("index {index} is out of bounds [0,2]")
        }
    }
}

impl Iterator for Vector {
    type Item = Length;

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

impl Vector {

    pub fn new(x: Length, y: Length, z: Length) -> Self { Self { x, y, z } }

    pub fn from<T>(x: f32, y: f32, z: f32) -> Self
    where
        T: uom::Conversion<f32, T = f32> + uom::si::length::Unit,
    {
        Self {
            x: Length::new::<T>(x),
            y: Length::new::<T>(y),
            z: Length::new::<T>(z),
        }
    }

    pub fn magnitude(&self) -> Length {
        let &Self { x, y, z } = self;
        (x*x + y*y + z*z).sqrt()
    }

    pub fn argmin(self) -> (usize, Length) { todo!() }
    pub fn norm(self) -> Length { todo!() }
    pub fn normalize(self) -> Self { todo!() }

}

#[cfg(test)]
mod tests {
    use crate::Vector;
    use float_eq::assert_float_eq;
    const EPS: f32 = f32::EPSILON;
    use uom::si::length::meter;
    use crate::uom::{nm, mm, cm, assert_uom_eq};
    use rstest::rstest;

    #[test]
    fn vector_components() {
        let v = Vector::new(mm(10.0), nm(1000.0), mm(2.0));
        assert_eq!(           v.x, mm(10.0  ));
        assert_uom_eq!(meter, v.y, mm( 0.001), r2nd <= EPS);
        assert_uom_eq!(meter, v.z, cm( 0.2  ), r2nd <= EPS);
    }

    #[test]
    fn mul_f32_for_vector() {
        let v = Vector::new(mm(1.0), mm(2.0), mm(3.0));
        let e = Vector::new(cm(1.0), cm(2.0), cm(3.0));
        let r = v * 10.0;
        assert_uom_eq!(meter, r.x, e.x, ulps <= 2);
        assert_uom_eq!(meter, r.y, e.y, ulps <= 2);
        assert_uom_eq!(meter, r.z, e.z, ulps <= 2);
    }

    #[rstest(/**/ x,  y,  z,  magnitude,
             case(0.0,  0.0,  0.0,  0.0),
             case(1.0,  0.0,  0.0,  1.0),
             case(0.0,  1.0,  0.0,  1.0),
             case(0.0,  0.0,  1.0,  1.0),
             case(3.0,  4.0,  0.0,  5.0),
             case(0.0, -3.0,  4.0,  5.0),
             case(5.0,  0.0, 12.0, 13.0),
             case(3.0,  4.0,  5.0,  7.0710678),
    )]
    fn vector_magnitude(x: f32, y: f32, z: f32, magnitude: f32) {
        let v = Vector::from::<meter>(x, y, z);
        assert_eq!(v.magnitude().get::<meter>(), magnitude);
    }

    #[test]
    fn index_for_vector() {
        let v = Vector::from::<meter>(1.0, 2.0, 3.0);
        assert_uom_eq!(meter, v[0], cm(100.0), ulps <= 1);
        assert_uom_eq!(meter, v[1], cm(200.0), ulps <= 1);
        assert_uom_eq!(meter, v[2], cm(300.0), ulps <= 1);
    }

}
