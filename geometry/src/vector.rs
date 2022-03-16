use std::ops::{Index, Mul};
use uom::si::f32::Length;

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

#[cfg(test)]
mod tests {
    use crate::Vector;
    use float_eq::assert_float_eq;
    const EPS: f32 = f32::EPSILON;
    use uom::si::length::meter;
    use crate::uom::{nm, mm, cm, assert_uom_eq};

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
        assert_uom_eq!(meter, r.x, e.x, ulps <= 1);
        assert_uom_eq!(meter, r.y, e.y, ulps <= 1);
        assert_uom_eq!(meter, r.z, e.z, ulps <= 1);
    }
}
