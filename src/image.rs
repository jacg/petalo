use crate::types::{Intensityf32, Index1_u, Index3_u};
use crate::weights::FOV;

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

// --------------------------------------------------------------------------------
//                  Conversion between 1d and 3d indices

use std::ops::{Add, Div, Mul, Rem};

pub fn index3_to_1<T>([ix, iy, iz]: [T; 3], [nx, ny, _nz]: [T; 3]) -> T
where
    T: Mul<Output = T> + Add<Output = T>
{
    ix + (iy + iz * ny) * nx
}

#[allow(clippy::many_single_char_names)]
pub fn index1_to_3<T>(i: T, [nx, ny, _nz]: [T; 3]) -> [T; 3]
where
    T: Mul<Output = T> +
    Div<Output = T> +
    Rem<Output = T> +
    Copy
{
    let z = i / (nx * ny);
    let r = i % (nx * ny);
    let y = r / nx;
    let x = r % nx;
    [x,y,z]
}


#[cfg(test)]
mod test_index_conversion {
    use super::*;
    use rstest::rstest;
    use crate::types::Index3_u;

    // -------------------- Some hand-picked examples ------------------------------
    #[rstest(/**/    size   , index3 , index1,
             // 1-d examples
             case([ 1, 1, 1], [0,0,0],   0),
             case([ 9, 1, 1], [3,0,0],   3),
             case([ 1, 8, 1], [0,4,0],   4),
             case([ 1, 1, 7], [0,0,5],   5),
             // Counting in binary: note digit reversal
             case([ 2, 2, 2], [0,0,0],   0),
             case([ 2, 2, 2], [1,0,0],   1),
             case([ 2, 2, 2], [0,1,0],   2),
             case([ 2, 2, 2], [1,1,0],   3),
             case([ 2, 2, 2], [0,0,1],   4),
             case([ 2, 2, 2], [1,0,1],   5),
             case([ 2, 2, 2], [0,1,1],   6),
             case([ 2, 2, 2], [1,1,1],   7),
             // Relation to decimal: note reversal
             case([10,10,10], [1,2,3], 321),
             case([10,10,10], [7,9,6], 697),
    )]
    fn hand_picked(size: Index3_u, index3: Index3_u, index1: usize) {
        assert_eq!(index3_to_1(index3, size), index1);
        assert_eq!(index1_to_3(index1, size), index3);
    }

    // -------------------- Exhaustive roundtrip testing ------------------------------
    use proptest::prelude::*;

    // A strategy that picks 3-d index limits, and a 1-d index guaranteed to lie
    // within those bounds.
    fn size_and_in_range_index() -> impl Strategy<Value = (Index3_u, usize)> {
        [1..200_usize, 1..200_usize, 1..200_usize]
            .prop_flat_map(|i| (Just(i), 1..(i[0] * i[1] * i[2])))
    }

    proptest! {
        #[test]
        fn index_roundtrip((size, index) in size_and_in_range_index()) {
            let there = index1_to_3(index, size);
            let back  = index3_to_1(there, size);
            assert_eq!(back, index)
        }

    }
}
