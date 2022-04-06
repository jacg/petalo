use crate::types::Weightf32;

#[allow(non_camel_case_types)] pub type Index1_u = usize;
#[allow(non_camel_case_types)] pub type Index3_u = [usize; 3];
#[allow(non_camel_case_types)] pub type BoxDim_u = [usize; 3];
pub type BoxDim = [LengthU; 3];

pub type Index3Weightf32 = (Index3_u, Weightf32);

pub type LengthI = geometry::uom::uomcrate::si::i32  ::Length;
pub type LengthU = geometry::uom::uomcrate::si::usize::Length;


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
