///! Test for geometry being used in conjunction with uom.


//use uom::fmt::DisplayStyle::Abbreviation;
use uom::si::f32::{Length, Time, Velocity};
use uom::si::{length  ::{nanometer, millimeter, centimeter},
              time    ::{nanosecond, picosecond},
              velocity:: meter_per_second};

// Making values from float literals seems to be very long-winded, so provide
// some pithily-named convenience constructors. These would probably have to be
// packed up in a constructor module in real life.
pub fn cm (x: f32) -> Length   {  Length::new::      <centimeter>(x) }
pub fn mm (x: f32) -> Length   {  Length::new::      <millimeter>(x) }
pub fn nm (x: f32) -> Length   {  Length::new::      < nanometer>(x) }
pub fn ns (x: f32) -> Time     {    Time::new::      <nanosecond>(x) }
pub fn ps (x: f32) -> Time     {    Time::new::      <picosecond>(x) }
pub fn m_s(x: f32) -> Velocity {Velocity::new::<meter_per_second>(x) }




#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Point, Vector};
    use float_eq::assert_float_eq;
    const EPS: f32 = f32::EPSILON;
    //use pretty_assertions::assert_eq;

    use uom::si::length::meter as m;

    macro_rules! assert_uom_eq {
        ($unit:ident, $lhs:expr, $rhs:expr, $algo:ident <= $tol:expr) => {
            assert_float_eq!($lhs.get::<$unit>(), $rhs.get::<$unit>(), $algo <= $tol)
        };
    }

    #[test]
    fn point_components() {
        let p = Point::<Length>::new(mm(10.0), nm(1000.0), mm(2.0));
        assert_eq!(       p.x, mm(10.0  ));
        assert_uom_eq!(m, p.y, mm( 0.001), r2nd <= EPS);
        assert_uom_eq!(m, p.z, cm( 0.2  ), r2nd <= EPS);
    }

    #[test]
    fn point_minus_point_same_unit() {
        // Difference between Points is a Vector
        let lhs      = Point ::<Length>::new(cm(3.0), mm( 20.0), cm( 8.0));
        let rhs      = Point ::<Length>::new(cm(2.0), cm(  4.0), mm(20.0));
        let expected = Vector::<Length>::new(cm(1.0), mm(-20.0), mm(60.0));
        let result: Vector<Length> = lhs - rhs;

        assert_uom_eq!(millimeter, result.x, expected.x, ulps <= 1);
        assert_uom_eq!(millimeter, result.y, expected.y, ulps <= 2);
        assert_uom_eq!(millimeter, result.z, expected.z, ulps <= 2);

    }

}

