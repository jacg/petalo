use units::todo::Intensityf32;

use crate::{
    fov::FOV,
    Index1_u, Index3_u,
    index::index3_to_1,
};

pub type ImageData = Vec<Intensityf32>;


#[derive(Clone)]
pub struct Image {
    pub fov: FOV,
    pub data: ImageData,
}



impl core::ops::IndexMut<Index1_u> for Image {
    #[inline]
    fn index_mut(&mut self, i: Index1_u) -> &mut Self::Output { &mut self.data[i] }
}

impl core::ops::Index<Index1_u> for Image {
    type Output = Intensityf32;
    #[inline]
    fn index(&self, i: Index1_u) -> &Self::Output { &self.data[i] }
}

impl core::ops::IndexMut<Index3_u> for Image {
    fn index_mut(&mut self, i3: Index3_u) -> &mut Self::Output {
        let i1 = index3_to_1(i3, self.fov.n);
        &mut self.data[i1]
    }
}

impl core::ops::Index<Index3_u> for Image {
    type Output = Intensityf32;
    fn index(&self, i3: Index3_u) -> &Self::Output {
        let i1 = index3_to_1(i3, self.fov.n);
        &self.data[i1]
    }
}

