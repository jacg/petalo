// Make uom crate reachable from our root crate, without clashing with our
// geometry::uom. TODO: this is a hack.
pub use uom as uomcrate;

//use uom::fmt::DisplayStyle::Abbreviation;
pub use uom::si::Quantity;
pub use uom::si::f32::{Length, Time, Velocity, Ratio};
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

pub fn ratio(x: f32) -> Ratio  {   Ratio::new::<uom::si::ratio::ratio>(x) }

// Reverse direction of the above. Rethink nomenclature once the dust has
// settled after the transition to uom is complete.
pub fn mm_(x: Length) -> f32 { x.get::<millimeter>() }

#[macro_export]
macro_rules! in_base_unit {
  ($value:expr) => {
    Quantity {
      dimension: std::marker::PhantomData,
      units: std::marker::PhantomData,
      value: $value,
    }
  };
}


#[cfg(test)]
use float_eq::assert_float_eq;

#[cfg(test)]
macro_rules! assert_uom_eq {
  ($unit:ident, $lhs:expr, $rhs:expr, $algo:ident <= $tol:expr) => {
    assert_float_eq!($lhs.get::<$unit>(), $rhs.get::<$unit>(), $algo <= $tol)
  };
}

#[cfg(test)]
pub (crate) use assert_uom_eq;


#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_name() {
    let v = vec![mm(1.0), cm(1.0)];
    let total: Length = v.into_iter().sum();
    assert_uom_eq!(nanometer, total, mm(11.0), ulps <= 1);
  }
}
