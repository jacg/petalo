/// Units which are simply type aliases for `f32` rather than having an
/// implementation as a `uom` `Quantity`.
///
/// This may be because:
///
/// + We do not know how to implement them in `uom`.
///
/// + There are other complications in the client code which make using `uom`
///   difficult, so we use plain `f32`s, but still want some clues in the source
///   as to what they represent.
///
/// + They might not be difficult to implement, but we haven't had the time to
///   think about it.

pub type Lengthf32    = f32;
pub type Timef32      = f32;
pub type Weightf32    = f32; // TODO uom Weight
pub type Ratiof32     = f32;
pub type Energyf32    = f32; // TODO uom Energy
pub type Chargef32    = f32; // TODO uom Charge
pub type Intensityf32 = f32; // TODO uom Intensity
