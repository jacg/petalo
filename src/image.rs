use std::path::Path;
use units::todo::Intensityf32;

use crate::{
    fov::FOV,
    Index1_u, Index3_u,
    index::index3_to_1,
    io,
};

pub type ImageData = Vec<Intensityf32>;


#[derive(Clone)]
pub struct Image {
    pub fov: FOV,
    pub data: ImageData,
}

impl Image {

    pub fn from_raw_file(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        Ok((&crate::io::raw::Image3D::read_from_file(path)?).into())
    }

    pub fn write_to_raw_file(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        io::raw::Image3D::from(self).write_to_file(path)?;
        Ok(())
    }

    pub fn ones(fov: FOV) -> Self {
        let [x,y,z] = fov.n;
        let size = x * y * z;
        Self { data: vec![1.0; size], fov}
    }

    pub fn new(fov: FOV, data: ImageData) -> Self {
        let [x, y, z] = fov.n;
        if data.len() != x * y * z {
            // TODO change panic to Option or Result
            panic!("Image data does not match dimensions {:?}", fov.n);
        };
        Image { fov, data }
    }

    pub fn empty(fov: FOV) -> Self {
        let [x,y,z] = fov.n;
        Self::new(fov, vec![0.0; x*y*z])
    }

    pub fn inverted(&self) -> Self {
        let mut inverted = self.clone();
        for e in inverted.data.iter_mut() { *e = 1.0 / *e }
        inverted
    }

    // A new empty data store with matching size
    pub fn zeros_buffer(fov: FOV) -> ImageData { let [x,y,z] = fov.n; vec![0.0; x*y*z] }

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

