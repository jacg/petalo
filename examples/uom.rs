/// UOM (Units Of Measure) is a crate which distinguishes different physical
/// units at type-level. For example, an f32 representing length, is not the
/// same type as an f32 representing time.
///
/// You can run this example with
/// ```
/// cargo run --example uom
/// ```
/// By default, the code which will generate compile-time errors is disabled. To
/// see these errors:
/// ```
/// cargo run --example uom --features compile-error
/// ```

use uom::fmt::DisplayStyle::Abbreviation;
use uom::si::f32::{Length, Time, Velocity};
use uom::si::{length  ::{kilometer, meter, millimeter},
              time    :: second,
              velocity:: meter_per_second};

// Making values from float literals seems to be very long-winded, so provide
// some pithily-named convenience constructors. These would probably have to be
// packed up in a constructor module in real life.
fn km(x: f32) -> Length { Length::new::< kilometer> (x) }
fn  m(x: f32) -> Length { Length::new::<     meter> (x) }
fn mm(x: f32) -> Length { Length::new::<millimeter> (x) }
fn  s(x: f32) -> Time   {   Time::new::<second>(x) }


fn main() {
    let c: Velocity = Velocity::new::<meter_per_second>(299_792_458.0);

    let l1: Length = km (3.2);         // Explicit type: Length
    let l2         =  m(40.2);         // Implicit type: Length
    let l3         = mm(99.2);
    let t          =  s(13.3);
    let v1: Velocity = (l1 + l2) / t;  // Length / Time -> Velocity
    let v2           =       l2  / t;

    #[cfg(feature = "compile-error")]
    let type_error = l1 + t;           // Cannot add length to time

    let two: f32 = 2.0;
    let two_l1 = two * l1;             // Multiplication by dimensionless number

    // Reusable format arguments
    let fmm = Length::format_args(millimeter, Abbreviation);
    let  fm = Length::format_args(     meter, Abbreviation);
    let fkm = Length::format_args( kilometer, Abbreviation);

    // Print out values along with their units
    println!("{} = {} = {} / 2",  fm.with(l1), fkm.with(l1), fm.with(two_l1));
    println!("{} = {}", fmm.with(l3),  fm.with(l3));
    println!("{}", v1.into_format_args(meter_per_second, Abbreviation));
    println!("{:?}", v2 / c);

    // Other example sets, defined below
    sep();    quantities();
    sep();    ncollide();
    sep();    parry();
}

// ===========================================================================
// Example showing how to create a set of `Quantity` type aliases for a
// different set of base units.

use uom::si::length::centimeter;

mod cgs {
    use uom::{ISQ, system};
    ISQ!(uom::si, f32, (millimeter, gram, microsecond, ampere, kelvin, mole, candela));
}

fn quantities() {
    let lm = uom::si::f32::Length::new::     <meter>(1.2345);
    let lc = cgs::         Length::new::<centimeter>(123.45);
    let t  = uom::si::f32::Time  ::new::    <second>(0.0015);

    println!("lm {}: {:?}"  , uom::si::length::description(), lm);
    println!("lc {}: {:?}\n", uom::si::length::description(), lc);

    println!("{:?} + {:?} = {:?}"  , lm, lc, (lm + lc));
    println!("{:?} + {:?} = {:?}\n", lc, lm, (lc + lm));

    println!("{:?} - {:?} = {:?}"  , lm, lc, (lm - lc));
    println!("{:?} - {:?} = {:?}\n", lc, lm, (lc - lm));

    println!("{:?} / {:?} = {:?}"  , lm,  t, (lm /  t));
    println!("{:?} / {:?} = {:?}"  , lc,  t, (lc /  t));
}

// ============================================================================
// Using UOM units in ncollide vectors/points works up to a point. Here we show
// that arithmetic ops on Point and Vector work, but Vector::norm() does not.

fn ncollide() {

    use ncollide3d as nc;
    use nc::math::{Vector, Point};

    type PL = nc::math::Point<Length>;
    type Pf = nc::math::Point<f32>;

    let p1f: Point<f32>    = Pf::new(  1.2,    1.3,    1.4);
    let p2f                = Pf::new(  2.1,    3.1,    4.1);

    let p1u: Point<Length> = PL::new(m(1.2), m(1.3), m(1.4));
    let p2u                = PL::new(m(2.1), m(3.1), m(4.1));

    let vf: Vector<f32>    = p2f - p1f;
    let vu: Vector<Length> = p2u - p1u;

    println!("{:?}", p2u - p1u);
    println!("{:?}", p2f - p1f);

    println!("{:?}", vf);
    println!("{:?}", vu);

    println!("{:?}", vf + vf);
    println!("{:?}", vu + vu);

    #[cfg(feature = "compile-error")]
    println!("{:?}", vu.norm());
    println!("{:?}", vf.norm());
}
// ============================================================================
// Using UOM units in parry3d vectors/points works up to a point. Here we show
// that arithmetic ops on Point and Vector work, but Vector::norm() does not.

fn parry() {
    use parry3d as p3d;

    type PL = p3d::math::Point<Length>;
    type Pf = p3d::math::Point<f32>;

    use p3d::math::{Vector, Point};

    let p1f: Point<f32>    = Pf::new(  1.2,    1.3,    1.4);
    let p2f                = Pf::new(  2.1,    3.1,    4.1);

    let p1u: Point<Length> = PL::new(m(1.2), m(1.3), m(1.4));
    let p2u                = PL::new(m(2.1), m(3.1), m(4.1));

    let vf: Vector<f32>    = p2f - p1f;
    let vu: Vector<Length> = p2u - p1u;

    println!("{:?}", p2u - p1u);
    println!("{:?}", p2f - p1f);

    println!("{:?}", vf);
    println!("{:?}", vu);

    println!("{:?}", vf + vf);
    println!("{:?}", vu + vu);

    #[cfg(feature = "compile-error")]
    println!("{:?}", vu.norm());
    println!("{:?}", vf.norm());

}
// ============================================================================
fn sep() {
    println!("================================================================================");
}
