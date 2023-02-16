//! Find indices and weights of the voxels coupled to a single LOR
//!
//! The algorithm is centred around two key simplifications:
//!
//! 1. Express the voxel size in terms of the components of the LOR's direction
//!    vector. This allows trivial calculation of how far we must move along the
//!    LOR before reaching a voxel boundary, in any dimension.
//!
//! 2. Exploit symmetry to simplify dealing with directions: flip axes so that
//!    the direction of the LOR has non-negative components. The algorithm can
//!    then assume that all progress is in the positive direction. Any voxels
//!    indices calculated by the algorithm, must be flipped back to the original
//!    coordinate system.

use units::todo::Weightf32;

use crate::{
    Index1_u,
    image::{Image, ImageData},
};

pub mod siddon;
pub use siddon::Siddon;


// Data needed by project_one_lor, both as input and output, because of the
// constrains imposed by `fold`
pub struct FoldState<'img, T: ?Sized> {
    pub backprojection: ImageData,
    pub system_matrix_row: SystemMatrixRow,
    pub image: &'img Image,
    pub projector_data: T,
}

pub type Fs<'i, S> = FoldState<'i, S>;

// ---------------------- Implementation -----------------------------------------
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

// ----- A trait for implementations of (geometrical?) system matrix calculations --------------------
pub trait SystemMatrix {

    //type Data;

    fn update_system_matrix_row(
        system_matrix_row: &mut SystemMatrixRow,
        lor: &LOR,
        fov:  FOV,
        data: &Self, //::Data,
    );

    // Sparse storage of the slice through the system matrix which corresponds
    // to the current LOR. Allocating these anew for each LOR had a noticeable
    // runtime cost, so we create them up-front and reuse them.
    // This should probably have a default implementation
    fn buffers(fov: FOV) -> SystemMatrixRow;
}

use crate::{
    LOR,
    fov::FOV,
};
