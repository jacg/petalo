use units::{Area, Length, Ratio, mm};
use units::uom::typenum::P2;
use crate::{Point, RatioVec, Dot};

/// Compute the length of the intersection of the line passing through points
/// `p1` and `p2` with a cylinder of radius `r` whose axis coincides with the
/// z-axis.
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
    let four_a_c = 4. * a * c;
    if b_squared <= four_a_c { return mm(0.0) }
    //
    let delta_t: Ratio = (b_squared - four_a_c).sqrt() / a;
    (delta_t * v).norm()
}

fn squared(l: Length) -> Area { l.powi(P2::new()) }

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use units::mm_;
    use units::float_eq::assert_float_eq;

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
}
