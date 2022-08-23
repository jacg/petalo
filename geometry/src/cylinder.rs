use units::{Area, Length, Length4, Ratio, ratio};
use units::uom::typenum::P2;
use crate::{Point, RatioVec, Dot};

/// Compute the length of the intersection of the line passing through points
/// `p1` and `p2` with a cylinder of radius `r` whose axis coincides with the
/// z-axis.
#[cfg(test)]
pub fn cylinder_line_intersection_length(p1: Point, p2: Point, r: Length) -> Length {
    // TODO: explain implementation
    let z = RatioVec::new(0., 0., 1.);
    let v = p2 - p1;
    let w = p1 - Point::zero();
    // Viète coefficients
    let a: Area = v.dot(v) - squared(v.dot(z));
    let b: Area = 2. * (v.dot(w) - v.dot(z) * w.dot(z));
    let c: Area = w.dot(w) - squared(w.dot(z)) - squared(r);
    // Check discriminant to see if line missed cylinder
    let b_squared = b * b;
    (delta_t(a, b_squared, c) * v).norm()
}

fn delta_t(a: Area, b_squared: Length4, c: Area) -> Ratio {
    let four_a_c = 4. * a * c;
    if b_squared <= four_a_c { return ratio(0.0) }
    (b_squared - four_a_c).sqrt() / a
}

// An obvious implementation, but probably less efficient than it could be,
// because the two function calls it uses will repeat many of the same
// calculations.
#[cfg(test)]
fn hollow_cylinder_line_intersection_length_slow(p1: Point, p2: Point, r: Length, dr: Length) -> Length {
    let inner = cylinder_line_intersection_length(p1, p2, r);
    let outer = cylinder_line_intersection_length(p1, p2, r+dr);
    outer - inner
}

// A (hopefully) more efficient implementation of the above. TODO: benchmark
pub fn hollow_cylinder_line_intersection_length(p1: Point, p2: Point, r: Length, dr: Length) -> Length {
    let z = RatioVec::new(0., 0., 1.);
    let v = p2 - p1;
    let w = p1 - Point::zero();
    // Viète coefficients
    let a:  Area = v.dot(v) - squared(v.dot(z));
    let b:  Area = 2. * (v.dot(w) - v.dot(z) * w.dot(z));
    let ci: Area = w.dot(w) - squared(w.dot(z)) - squared(r   );
    let co: Area = w.dot(w) - squared(w.dot(z)) - squared(r+dr);

    let b_squared = b * b;

    let delta_ti = delta_t(a, b_squared, ci);
    let delta_to = delta_t(a, b_squared, co);
    let dt = delta_to - delta_ti;
    (dt * v).norm()
}

fn squared(l: Length) -> Area { l.powi(P2::new()) }

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use proptest::prelude::*;
    use units::uom::si::length::millimeter;
    use units::{mm, mm_};
    use units::float_eq::assert_float_eq;
    use units::assert_uom_eq;

    #[rstest(/**/   x1,    y1,    z1,     x2,   y2,    z2,     r, expected,
             case( 10.0, -10.0,  0.0,   10.0, 10.0,   0.0,   1.0, 0.0), // miss completely on right
             case(-10.0,  10.0,  0.0,   10.0, 10.0,   0.0,   2.0, 0.0), // miss completely above
             case( 10.0, -10.0,  9.0,   10.0, 10.0,  19.0,   1.0, 0.0), // as above, but ...
             case(-10.0,  10.0, -3.0,   10.0, 10.0,  19.0,   2.0, 0.0), // ... different zs
             case(  0.0, -10.0,  0.0,    0.0, 10.0,   0.0,   3.0, 6.0), // along vertical   diameter
             case(-10.0,   0.0,  0.0,   10.0,  0.0,   0.0,   4.0, 8.0), // along horizontal diameter
             case(-10.0,   0.0, 10.0,   10.0,  0.0, -10.0,   4.0, 8.0 * 2.0_f32.sqrt()), // at 45 degrees to z-axis
             case(-10.0,   2.5,  0.0,   10.0,  2.5,   0.0,   5.0, 5.0 * 3.0_f32.sqrt()), // off-z-axis by r/2
    )]
    fn test_cylinder_line_intersection_length(
        x1: f32, y1: f32, z1: f32,
        x2: f32, y2: f32, z2: f32,
        r: f32,
        expected: f32,
    ) {
        let p1 = Point::new(mm(x1), mm(y1), mm(z1));
        let p2 = Point::new(mm(x2), mm(y2), mm(z2));
        let calculated = mm_(cylinder_line_intersection_length(p1, p2, mm(r)));
        assert_float_eq!(calculated, expected, ulps <= 1);
    }

    proptest! {
        #[test]
        fn hollow_cylinder_intersection_length_compare_implementations(
            x1 in    0.1..(100.0 as f32),
            y1 in -100.0..(100.0 as f32),
            z1 in -100.0..(100.0 as f32),
            x2 in -100.0..(  0.0 as f32),
            y2 in -100.0..(100.0 as f32),
            z2 in -100.0..(100.0 as f32),
            r  in  100.0..(300.0 as f32),
            dr in    1.0..( 30.0 as f32),
        ) {
            let p1 = Point::new(mm(x1), mm(y1), mm(z1));
            let p2 = Point::new(mm(x2), mm(y2), mm(z2));
            let (r, dr) = (mm(r), mm(dr));
            let slow = hollow_cylinder_line_intersection_length_slow(p1, p2, r, dr);
            let fast = hollow_cylinder_line_intersection_length     (p1, p2, r, dr);
            assert_uom_eq!(millimeter, fast, slow, rel <= 0.0001);
        }
    }
}
