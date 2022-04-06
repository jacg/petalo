use crate::types::Weightf32;

#[allow(non_camel_case_types)] pub type Index1_u = usize;
#[allow(non_camel_case_types)] pub type Index3_u = [usize; 3];
#[allow(non_camel_case_types)] pub type BoxDim_u = [usize; 3];
pub type BoxDim = [LengthU; 3];

pub type Index3Weightf32 = (Index3_u, Weightf32);

pub type LengthI = geometry::uom::uomcrate::si::i32  ::Length;
pub type LengthU = geometry::uom::uomcrate::si::usize::Length;
