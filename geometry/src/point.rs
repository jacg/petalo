use std::ops::{AddAssign, Index, Sub, IndexMut};
use crate::uom::mmps::f32::Length;
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
    fn add_assign(&mut self, delta: Vector) {
        self.x += delta.x;
        self.y += delta.y;
        self.z += delta.z;
    }
}

impl Index<usize> for Point {
    type Output = Length;
    fn index(&self, index: usize) -> &Length {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("index {index} is out of bounds [0,2]")
        }
    }
}

impl IndexMut<usize> for Point {
    fn index_mut(&mut self, index: usize) -> &mut Length {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("index {index} is out of bounds [0,2]")
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Point, Vector};
    const EPS: f32 = f32::EPSILON;
    use uom::si::length::{millimeter, meter};
    use crate::uom::{nm, mm, cm, assert_uom_eq};

    #[test]
    fn point_components() {
        let p = Point::new(mm(10.0), nm(1000.0), mm(2.0));
        assert_eq!(           p.x, mm(10.0  ));
        assert_uom_eq!(meter, p.y, mm( 0.001), r2nd <= EPS);
        assert_uom_eq!(meter, p.z, cm( 0.2  ), r2nd <= EPS);
    }

    #[test]
    fn sub_for_point() {
        // Difference between Points is a Vector
        let lhs      = Point ::new(cm(3.0), mm( 20.0), cm( 8.0));
        let rhs      = Point ::new(cm(2.0), cm(  4.0), mm(20.0));
        let expected = Vector::new(cm(1.0), mm(-20.0), mm(60.0));
        let result: Vector = lhs - rhs;

        assert_uom_eq!(millimeter, result.x, expected.x, ulps <= 1);
        assert_uom_eq!(millimeter, result.y, expected.y, ulps <= 2);
        assert_uom_eq!(millimeter, result.z, expected.z, ulps <= 2);
    }

    #[test]
    fn sub_for_ref_point() {
        // Difference between Points is a Vector
        let lhs      = Point ::new(cm(3.0), mm( 20.0), cm( 8.0));
        let rhs      = Point ::new(cm(2.0), cm(  4.0), mm(20.0));
        let expected = Vector::new(cm(1.0), mm(-20.0), mm(60.0));
        let result: Vector = &lhs - &rhs;

        assert_uom_eq!(millimeter, result.x, expected.x, ulps <= 1);
        assert_uom_eq!(millimeter, result.y, expected.y, ulps <= 2);
        assert_uom_eq!(millimeter, result.z, expected.z, ulps <= 2);
    }

    #[test]
    fn addassign_for_point() {
        let mut p = Point ::new(cm( 1.0), cm( 2.0), cm(3.0));
        let     v = Vector::new(mm(10.0), mm(15.0), cm(2.5));
        let xpct  = Point ::new(cm( 2.0), mm(35.0), cm(5.5));
        p += v;
        assert_uom_eq!(meter, p.x, xpct.x, ulps <= 1);
        assert_uom_eq!(meter, p.y, xpct.y, ulps <= 1);
        assert_uom_eq!(meter, p.z, xpct.z, ulps <= 1);
    }

    #[test]
    fn index_for_point_in_bounds() {
        let p = Point ::new(cm(1.0), cm(2.0), cm(3.0));
        assert_eq!(p[0], p.x);
        assert_eq!(p[1], p.y);
        assert_eq!(p[2], p.z);
        assert_eq!(p[0], cm(1.0));
        assert_eq!(p[1], cm(2.0));
        assert_eq!(p[2], cm(3.0));
    }

    #[test]
    #[should_panic]
    fn index_for_point_out_of_bounds() {
        let p = Point ::new(cm(1.0), cm(2.0), cm(3.0));
        p[3];
    }

    #[test]
    fn index_mut_for_point_in_bounds() {
        let mut p = Point ::new(cm(1.0), cm(2.0), cm(3.0));
        p[0] = cm(4.0);
        p[1] = cm(5.0);
        p[2] = cm(6.0);
        assert_eq!(p.x, cm(4.0));
        assert_eq!(p.y, cm(5.0));
        assert_eq!(p.z, cm(6.0));
    }

    #[test]
    #[should_panic]
    fn index_mut_for_point_out_of_bounds() {
        let mut p = Point ::new(cm(1.0), cm(2.0), cm(3.0));
        p[3] = cm(4.0);
    }

}
