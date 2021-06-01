use crate::io::raw;
use crate::types::{Length, Point, Intensity, Ratio};
use crate::mlem::{Image, ImageData};
use crate::weights::VoxelBox;

type BoxErr<T> = Result<T, Box<dyn std::error::Error>>;

pub fn load_image(filename: &std::path::PathBuf, vbox: VoxelBox) -> BoxErr<Image> {
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


#[derive(Clone)]
pub enum ROI {
    Sphere((Length, Length, Length), Length),
    CylinderZ((Length, Length), Length),
}

// TODO replace vec with iterator in output
impl Image {

    pub fn values_inside_roi(&self, roi: ROI) -> ImageData {
        let mut out = vec![];
        let roi_contains: Box<dyn Fn(Point) -> bool> = match roi {
            ROI::Sphere((cx, cy, cz), radius) => Box::new(move |p: Point| {
                let (x,y,z) = (p.x - cx, p.y - cy, p.z - cz);
                x*x + y*y + z*z < radius * radius
            }),

            ROI::CylinderZ((cx, cy), radius) => Box::new(move |p: Point| {
                let (x, y) = (p.x - cx, p.y - cy);
                x*x + y*y < radius*radius
            })
        };
        for (index, value) in self.data.iter().copied().enumerate() {
            let p = self.vbox.voxel_centre1(index);
            if roi_contains(p) { out.push(value) }
        }
        out
    }

}

#[cfg(test)]
mod test_in_roi {
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
        let inside = image.values_inside_roi(ROI::Sphere(centre, r));
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
        let inside = image.values_inside_roi(ROI::CylinderZ(centre, r));
        println!("{} {}", inside.len(), expected_len);
        assert!(inside.len() == expected_len);
    }

}

// TODO stop reinventing this wheel
fn mean(data: &ImageData) -> Option<Intensity> {
    data.iter().cloned().reduce(|a, b| a+b).map(|s| s / data.len() as Intensity)
}

impl Image {
    pub fn sphere_crcs(
        &self,
        rois           : &[ROI],        roi_activities: &[Intensity],
        background_rois: &[ROI], background_activity  :   Intensity,
    )
        -> Vec<Ratio>
    {
        let background_measured = background_rois.iter().cloned()
            .map(|roi| mean(&self.values_inside_roi(roi)).unwrap())
            .reduce(|a,b| a+b)
            .unwrap();

        let mut results = vec![];
        for (roi, roi_activity) in rois.iter().cloned().zip(roi_activities) {
            let roi_measured = mean(&self.values_inside_roi(roi)).unwrap();
            results.push(((roi_measured / background_measured) - 1.0) /
                         ((roi_activity / background_activity) - 1.0)  )
        }
        results
    }
}

#[cfg(test)]
mod test_crc {
    use super::*;

}
