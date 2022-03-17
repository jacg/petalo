// Make uom crate reachable from our root crate, without clashing with our
// geometry::uom. TODO: this is a hack.
pub use uom as uomcrate;

//use uom::fmt::DisplayStyle::Abbreviation;
pub use uom::si::f32::{Length, Time, Velocity};
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


#[allow(unused_macros)]
macro_rules! assert_uom_eq {
  ($unit:ident, $lhs:expr, $rhs:expr, $algo:ident <= $tol:expr) => {
    assert_float_eq!($lhs.get::<$unit>(), $rhs.get::<$unit>(), $algo <= $tol)
  };
}

#[allow(unused_imports)]
pub (crate) use assert_uom_eq;
