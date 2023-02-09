use crate::{
    LOR,
    fov::FOV,
    gauss::Gaussian,
    image::{ImageData, Image},
    system_matrix::SystemMatrixRow,
};

pub use siddon::Siddon;


pub mod siddon;

pub trait Projector {
    fn project_one_lor<'i, 'g>(fold_state: FoldState<'i, 'g>, lor: &LOR) -> FoldState<'i, 'g>;

    // Sparse storage of the slice through the system matrix which corresponds
    // to the current LOR. Allocating these anew for each LOR had a noticeable
    // runtime cost, so we create them up-front and reuse them.
    // This should probably have a default implementation
    fn buffers(fov: FOV) -> SystemMatrixRow;
}


pub type FoldState<'i, 'g> = (ImageData, SystemMatrixRow, &'i Image, &'g Option<Gaussian>);
