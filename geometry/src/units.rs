use uom::si::Dimension;
pub type InvertDimension<D> = uom::si::ISQ<
    <<D as Dimension>::L  as uom::lib::ops::Neg>::Output,
    <<D as Dimension>::M  as uom::lib::ops::Neg>::Output,
    <<D as Dimension>::T  as uom::lib::ops::Neg>::Output,
    <<D as Dimension>::I  as uom::lib::ops::Neg>::Output,
    <<D as Dimension>::Th as uom::lib::ops::Neg>::Output,
    <<D as Dimension>::N  as uom::lib::ops::Neg>::Output,
    <<D as Dimension>::J  as uom::lib::ops::Neg>::Output>;

pub mod mmps {

  use uom::si::{
    length::millimeter,
    mass::kilogram,
    time::picosecond,
    electric_current::ampere,
    thermodynamic_temperature::kelvin,
    amount_of_substance::mole,
    luminous_intensity::candela,
  };

  // TODO: replace with system! macro, once it has been fixed in uom
  type Units = dyn uom::si::Units<
      f32,
    length                    = millimeter,
    mass                      = kilogram,
    time                      = picosecond,
    electric_current          = ampere,
    thermodynamic_temperature = kelvin,
    amount_of_substance       = mole,
    luminous_intensity        = candela>;

  pub mod f32 {
    use uom::{ISQ, system, si::Quantity};
    ISQ!(uom::si, f32, (millimeter, kilogram, picosecond, ampere, kelvin, mole, candela));

    use uom::typenum::{P2, N1, Z0};
    pub type PerLength   = Quantity<super::super::InvertDimension<uom::si::length::Dimension>, super::Units, f32>;
    pub type AreaPerMass = Quantity<uom::si::ISQ<P2, N1, Z0, Z0, Z0, Z0, Z0>                 , super::Units, f32>;


    /// The full circle constant (τ) Equal to 2π.
    pub const TWOPI: Angle = Angle {
        dimension: std::marker::PhantomData,
        units: std::marker::PhantomData,
        value: std::f32::consts::TAU,
    };
  }

  pub mod i32 {
    use uom::{ISQ, system};
    ISQ!(uom::si, i32, (millimeter, kilogram, picosecond, ampere, kelvin, mole, candela));
  }

  pub mod usize {
    use uom::{ISQ, system};
    ISQ!(uom::si, usize, (millimeter, kilogram, picosecond, ampere, kelvin, mole, candela));
  }

}

//use uom::fmt::DisplayStyle::Abbreviation;
pub use uom::si::Quantity;
pub use mmps::f32::{Angle, TWOPI, Length, Time, Velocity, Ratio, Mass, PerLength, AreaPerMass};
mod units {
  pub use uom::si::{length  ::{nanometer, millimeter, centimeter},
                    time    ::{nanosecond, picosecond},
                    mass    ::kilogram,
                    velocity::meter_per_second,
                    ratio   ::ratio,
                    angle   ::{radian, revolution},
  };
}

// Making uom quantities from float literals is very long-winded, so provide
// some pithily-named convenience constructors. These would probably have to be
// packed up in a constructor module in real life.

/// Generate a pair of functions for converting between f32 and uom quantities.
///
/// wrap!(WRAP_NAME UNWRAP_NAME QUANTITY UNIT);
///
/// The wrapping function is called WRAP_NAME and returns QUANTITY by
/// interpreting its argument as UNIT. The function UNWRAP_NAME is the inverse
/// of WRAP_NAME.
macro_rules! wrap {
  ($wrap_name:ident $unwrap_name:ident $quantity:ident $unit:ident ) => {
    pub fn   $wrap_name(x: f32) -> $quantity { $quantity::new::<units::$unit>(x) }
    pub fn $unwrap_name(x: $quantity) -> f32 {          x.get::<units::$unit>( ) }
  };
}

wrap!(cm     cm_     Length         centimeter);
wrap!(mm     mm_     Length         millimeter);
wrap!(nm     nm_     Length          nanometer);
wrap!(ns     ns_     Time           nanosecond);
wrap!(ps     ps_     Time           picosecond);
wrap!(m_s    m_s_    Velocity meter_per_second);
wrap!(kg     kg_     Mass             kilogram);
wrap!(ratio  ratio_  Ratio               ratio);
wrap!(radian radian_ Angle              radian);
wrap!(turn   turn_   Angle          revolution);

pub fn mm_ps (x: f32) -> Velocity { m_s (x / m_s(1.0).value) }
pub fn mm_ps_(x: Velocity) -> f32 { m_s_(x * m_s(1.0).value) }


#[macro_export]
macro_rules! in_base_unit {
  ($value:expr) => {
    $crate::Quantity {
      dimension: std::marker::PhantomData,
      units: std::marker::PhantomData,
      value: $value,
    }
  };
}


#[macro_export]
macro_rules! assert_uom_eq {
  ($unit:ident, $lhs:expr, $rhs:expr, $algo:ident <= $tol:expr) => {
    float_eq::assert_float_eq!($lhs.get::<$unit>(), $rhs.get::<$unit>(), $algo <= $tol)
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
    use units::nanometer;
    assert_uom_eq!(nanometer, total, mm(11.0), ulps <= 1);
  }
}
