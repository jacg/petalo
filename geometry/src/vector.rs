use std::ops::{Index, IndexMut, Mul, Sub};
use units::{Area, Length, Ratio};

use units::{mm, ratio};

type NcVector = ncollide3d::math::Vector::<f32>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vec<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T: Sub<Output = T>> Sub for Vec<T> {
    type Output = Vec<<T as Sub>::Output>;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T, RHS> Mul<RHS> for Vec<T>
where
    T: Mul<RHS, Output = T>,
    RHS: Copy,
{
    type Output = Vec<<T as Mul<RHS>>::Output>;

    fn mul(self, rhs: RHS) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

pub type   Vector = Vec<Length>;
pub type RatioVec = Vec<Ratio>;
pub type  AreaVec = Vec<Area>;

impl Mul<Vector> for Ratio {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Self::Output {
        Vector {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl Mul<RatioVec> for Vector {
    type Output = Self;
    fn mul(self, rhs: RatioVec) -> Self::Output {
        Vector {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl Mul<Vector> for RatioVec {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Self::Output {
        Vector {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
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

impl IndexMut<usize> for Vector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
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

    pub fn zero() -> Self { Self::new(mm(0.), mm(0.), mm(0.)) }

    pub fn xyz<T>(x: f32, y: f32, z: f32) -> Self
    where
        T: units::uom::Conversion<f32, T = f32> + units::uom::si::length::Unit,
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

    #[allow(clippy::collapsible_else_if)]
    pub fn argmin(self) -> (usize, Length) {
        if self.x <= self.y {
            if self.x <= self.z { (0, self.x) }
            else                { (2, self.z) }
        } else {
            if self.y <= self.z { (1, self.y) }
            else                { (2, self.z) }
        }
    }

    pub fn norm        (self) -> Length { mm(NcVector::from(self).norm        ())           }
    pub fn norm_squared(self) -> Area   { mm(NcVector::from(self).norm_squared()) * mm(1.0) }

    pub fn normalize(self) -> RatioVec {
        let n = NcVector::from(self).normalize();
        let (x, y, z) = (ratio(n.x), ratio(n.y), ratio(n.z));
        RatioVec {x, y, z}
    }

    pub fn component_div(&self, rhs: RatioVec) -> Self {
        Vector {
            x: self.x / rhs.x,
            y: self.y / rhs.y,
            z: self.z / rhs.z,
        }
    }

    pub fn component_mul(&self, rhs: RatioVec) -> Self {
        Vector {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }

}


impl RatioVec{
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            x: ratio(x),
            y: ratio(y),
            z: ratio(z)
        }
    }
}



pub trait Dot<RHS = Self> {
    type Output;
    fn dot(&self, other: RHS) -> Self::Output;
}

impl Dot for Vector {
    type Output = Area;
    fn dot(&self, other: Vector) -> Area {
        self.x * other.x +
        self.y * other.y +
        self.z * other.z
    }
}

impl Dot<RatioVec> for Vector {
    type Output = Length;
    fn dot(&self, other: RatioVec) -> Length {
        self.x * other.x +
        self.y * other.y +
        self.z * other.z
    }
}

impl Dot<Vector> for RatioVec {
    type Output = Length;
    fn dot(&self, other: Vector) -> Length {
        self.x * other.x +
        self.y * other.y +
        self.z * other.z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const EPS: f32 = f32::EPSILON;
    use units::uom::si::area::square_millimeter;
    use units::uom::si::length::{meter, millimeter};
    use units::{nm, mm, cm, assert_uom_eq};
    use rstest::rstest;
    use proptest::prelude::*;

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

    #[allow(clippy::excessive_precision)]
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
        let v = Vector::xyz::<meter>(x, y, z);
        assert_eq!(v.magnitude().get::<meter>(), magnitude);
    }

    #[test]
    fn index_for_vector() {
        let v = Vector::xyz::<meter>(1.0, 2.0, 3.0);
        assert_uom_eq!(meter, v[0], cm(100.0), ulps <= 1);
        assert_uom_eq!(meter, v[1], cm(200.0), ulps <= 1);
        assert_uom_eq!(meter, v[2], cm(300.0), ulps <= 1);
    }

    proptest! {
        #[test]
        fn norm_equals_magnitude(
            x in -100.0..100.0_f32,
            y in -100.0..100.0_f32,
            z in -100.0..100.0_f32,
        ) {
            let v = Vector::new(mm(x), mm(y), mm(z));
            let n: Length = v.norm();
            let m: Length = v.magnitude();
            assert_uom_eq!(millimeter, n, m, ulps <=1);
        }
    }

    proptest! {
        #[test]
        fn norm_squared(
            x in -100.0..100.0_f32,
            y in -100.0..100.0_f32,
            z in -100.0..100.0_f32,
        ) {
            let v = Vector::new(mm(x), mm(y), mm(z));
            let n: Length = v.norm();
            let s: Area   = v.norm_squared();
            assert_uom_eq!(square_millimeter, n*n, s, ulps <=2);
        }
    }
}
