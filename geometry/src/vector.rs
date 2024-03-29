use std::ops::{Index, IndexMut, Add, Mul, Sub};
use units::{Area, Length, Ratio, in_base_unit};

use units::ratio;

type NcVector = ncollide3d::math::Vector::<f32>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vect<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Vect<T> {
    pub fn new(x: T, y: T, z: T) -> Self { Self { x, y, z } }
}

impl<T: units::uom::num_traits::Zero + Copy> Vect<T> {
    pub fn zero() -> Self { let zero = T::zero(); Self::new(zero, zero, zero) }
}

impl<LHS> Vect<LHS> {
    pub fn dot<RHS, Out>(self, other: Vect<RHS>) -> Out
    where
        LHS: Mul<RHS, Output=Out> + Copy,
        Out: Add<Output=Out>,
    {
        self.x * other.x +
        self.y * other.y +
        self.z * other.z
    }

    pub fn elementwise_mul<RHS, Out>(self, other: Vect<RHS>) -> Vect<Out>
    where
        LHS: Mul<RHS, Output=Out> + Copy,
    {
        Vect::<Out> {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }
    }
}

impl<T: Sub<Output = T>> Sub for Vect<T> {
    type Output = Vect<<T as Sub>::Output>;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T, RHS, Out> Mul<RHS> for Vect<T>
where
    T: Mul<RHS, Output=Out>,
    RHS: Copy,
{
    type Output = Vect<Out>;

    fn mul(self, rhs: RHS) -> Self::Output {
        Vect::<Out> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

pub type   Vector = Vect<Length>;
pub type RatioVec = Vect<Ratio>;
pub type _AreaVec = Vect<Area>;

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

impl<T> Index<usize> for Vect<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("index {index} is out of bounds [0,2]")
        }
    }
}

impl<T> IndexMut<usize> for Vect<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("index {index} is out of bounds [0,2]")
        }
    }
}

impl<T> Iterator for Vect<T> {
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        todo!()
    }
}

impl Vector {

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

    pub fn norm        (self) -> Length { in_base_unit!(NcVector::from(self).norm        ()) }
    pub fn norm_squared(self) -> Area   { in_base_unit!(NcVector::from(self).norm_squared()) }

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
