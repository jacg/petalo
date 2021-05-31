use crate::io::raw;
use crate::types::Length;
use crate::mlem::{Image, ImageData};
use crate::weights::VoxelBox;

type BoxErr<T> = Result<T, Box<dyn std::error::Error>>;

fn load_image(filename: &std::path::PathBuf, vbox: VoxelBox) -> BoxErr<Image> {
    let data = raw::read(filename)?.collect::<Result<_,_>>()?;
    Ok(Image::new(vbox, data)) // TODO: Upgrade Image::new from panic to Result
}

#[cfg(test)]
fn get_sample_image() -> BoxErr<Image> {
    let filename = std::path::PathBuf::from;
    let vbox = VoxelBox::new((180.0, 180.0, 180.0), (60, 60, 60));
    load_image(&filename("data/out/mlem/60_60_60_tof_OFF_00.raw"), vbox)
}

#[cfg(test)]
mod test_get_sample_image {
    use super::*;
    #[test]
    fn load_an_image_file() -> BoxErr<()> {
        let image = get_sample_image()?;
        assert!(image.data.len() == 60 * 60 * 60);
        Ok(())
    }
}



type Sphere = ((Length, Length, Length), Length);
type Cylinder = ((Length, Length), Length);

// TODO replace vec with iterator in output
impl Image {

    pub fn values_inside_sphere(&self, sphere: Sphere) -> ImageData {
        let mut out = vec![];
        let ((cx, cy, cz), radius) = sphere;
        for (index, value) in self.data.iter().copied().enumerate() {
            let p = self.vbox.voxel_centre1(index);
            let (x, y, z) = (p.x - cx, p.y - cy, p.z - cz);
            let d_squared = x*x + y*y + z*z;
            if d_squared < radius*radius { out.push(value); }
        }
        out
    }

    pub fn values_inside_cylinder(&self, cylinder: Cylinder) -> ImageData {
        let mut out = vec![];
        let ((cx, cy), radius) = cylinder;
        for (index, value) in self.data.iter().copied().enumerate() {
            let p = self.vbox.voxel_centre1(index);
            let (x, y) = (p.x - cx, p.y - cy);
            let d_squared = x*x + y*y;
            if d_squared < radius*radius { out.push(value); }
        }
        out
    }

}

#[cfg(test)]
mod test_in_sphere {
    use super::*;
    use rstest::rstest;

    // TODO these tests only verify the number of included voxels but do not
    // check whether the correct ones have been included

    // Arrange for outer voxels to be centred at +/- 100.0, when n = 10
    const MAGIC: Length = 10.0 / 9.0 * 200.0;

    #[rstest(/**/    l ,  n,        centre        ,    r , expected_len,
             case(MAGIC, 10, (  0.0,   0.0,   0.0), 173.3, 1000), // r > sqrt(3) * 100; all voxel centres inside sphere
             case(MAGIC, 10, (  0.0,   0.0,   0.0), 173.2,  992), // r < sqrt(3) * 100; 992 = 8 corners missing
             case(MAGIC,  9, (  0.0,   0.0,   0.0), 173.2,  729), // all 9^3 included: coarser grid => outer centres closer
             case(MAGIC, 10, (200.0, 200.0, 200.0), 173.3,    1), // single voxel at one corner of box
             case(MAGIC, 10, (200.0, 200.0, 200.0), 173.2,    0), // slightly smaller r excludes the corner
    )]
    fn number_of_included_voxels_in_sphere(
        l: Length,
        n: usize,
        centre: (Length, Length, Length),
        r: Length,
        expected_len: usize
    ) {
        let data = vec![1.0; n*n*n];
        let vbox = VoxelBox::new((l,l,l), (n,n,n));
        let image = Image::new(vbox, data);
        let inside = image.values_inside_sphere((centre, r));
        println!("{} {}", inside.len(), expected_len);
        assert!(inside.len() == expected_len);
    }

    #[rstest(/**/    l ,  n,        centre        ,    r , expected_len,
             case(MAGIC, 10, (  0.0,   0.0), 141.5, 1000), // r > sqrt(2) * 100; all voxel centres inside sphere
             case(MAGIC, 10, (  0.0,   0.0), 141.4,  960), // r < sqrt(2) * 100; 960 = 4*n corners missing
             case(MAGIC,  9, (  0.0,   0.0), 141.4,  729), // all 9^3 included: coarser grid => outer centres closer
             case(MAGIC, 10, (200.0, 200.0), 141.5,   10), // one row of voxels at one corner of box
             case(MAGIC, 10, (200.0, 200.0), 141.4,    0), // slightly smaller r excludes the corner
    )]
    fn number_of_included_voxels_in_z_cylinder(
        l: Length,
        n: usize,
        centre: (Length, Length),
        r: Length,
        expected_len: usize
    ) {
        let data = vec![1.0; n*n*n];
        let vbox = VoxelBox::new((l,l,l), (n,n,n));
        let image = Image::new(vbox, data);
        let inside = image.values_inside_cylinder((centre, r));
        println!("{} {}", inside.len(), expected_len);
        assert!(inside.len() == expected_len);
    }

}
