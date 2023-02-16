//! Calculation of system matrix elements for use in forward and backward
//! projections, both in MLEM iterations and sensitivity matrix constructions.

// ----- The trait --------------------------------------------------------------------

/// Interface for calculation of (geometrical components of?) system matrix elements
pub trait SystemMatrix {

    type Data;

    fn data(&self) -> Self::Data;

    /// Calculate the probabilities of a decay occurring in the voxels of `fov`
    /// being detected in `lor`. Place the results in the output parameter
    /// `system_matrix_row`.
    fn update_system_matrix_row(
        system_matrix_row: &mut SystemMatrixRow,
        lor: &LOR,
        fov:  FOV,
        data: &Self::Data,
    );

    // Sparse storage of the slice through the system matrix which corresponds
    // to the current LOR. Allocating these anew for each LOR had a noticeable
    // runtime cost, so we create them up-front and reuse them.
    // This should probably have a default implementation
    fn buffers(fov: FOV) -> SystemMatrixRow;
}

// ----- Implementations of the trait -----------------------------------------------
pub mod siddon;
pub use siddon::Siddon;

// ----- Storage of system matrix elements. Only one row is relevant at any single time ------
pub type SystemMatrixElement = (Index1_u, Weightf32);

pub struct SystemMatrixRow(pub Vec<SystemMatrixElement>);

impl SystemMatrixRow {
    pub fn iter(&self) -> std::slice::Iter<SystemMatrixElement> { self.0.iter() }
    pub fn clear(&mut self) { self.0.clear(); }
}

impl IntoIterator for SystemMatrixRow {
    type Item = SystemMatrixElement;
    type IntoIter = std::vec::IntoIter<Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a SystemMatrixRow {
    type Item = SystemMatrixElement;
    type IntoIter = std::iter::Cloned<std::slice::Iter<'a, Self::Item>>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter().cloned()
    }
}

// ----- Imports ------------------------------------------------------------------------------------------
use units::todo::Weightf32;

use crate::{
    LOR, Index1_u,
    fov::FOV,
};
