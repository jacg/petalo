use crate::io::raw;
use crate::types::{Length, Point, Intensity, Ratio};
use crate::mlem::{Image, ImageData};
use crate::weights::VoxelBox;

type BoxErr<T> = Result<T, Box<dyn std::error::Error>>;

pub fn load_image(filename: &std::path::Path, vbox: VoxelBox) -> BoxErr<Image> {
    let data = raw::read(filename)?.collect::<Result<_,_>>()?;
    Ok(Image::new(vbox, data)) // TODO: Upgrade Image::new from panic to Result
}


#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum ROI {
    Sphere((Length, Length, Length), Length),
    CylinderX((Length, Length), Length),
    CylinderY((Length, Length), Length),
    CylinderZ((Length, Length), Length),
}

impl ROI {
    pub fn contains_fn(self) -> Box<dyn Fn(Point) -> bool> {
        match self {
            ROI::Sphere((cx, cy, cz), radius) => Box::new(move |p: Point| {
                let (x,y,z) = (p.x - cx, p.y - cy, p.z - cz);
                x*x + y*y + z*z < radius * radius
            }),

            ROI::CylinderX((cy, cz), radius) => Box::new(move |p: Point| {
                let (y, z) = (p.y - cy, p.z - cz);
                y*y + z*z < radius*radius
            }),

            ROI::CylinderY((cx, cz), radius) => Box::new(move |p: Point| {
                let (x, z) = (p.x - cx, p.z - cz);
                x*x + z*z < radius*radius
            }),

            ROI::CylinderZ((cx, cy), radius) => Box::new(move |p: Point| {
                let (x, y) = (p.x - cx, p.y - cy);
                x*x + y*y < radius*radius
            }),

        }
    }
}

// TODO replace vec with iterator in output
impl Image {

    pub fn values_inside_roi(&self, roi: ROI) -> ImageData {
        let mut out = vec![];
        let roi_contains = roi.contains_fn();
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

    const I: Intensity =  1.0;
    const O: Intensity = -1.0;
    const R: Length = 1.0;
    use ROI::*;
    #[rstest(/**/  roi              , expected,
             // Should pick exactly one voxel
             case(Sphere((O,O,O), R),   1),
             case(Sphere((I,O,O), R),   2),
             case(Sphere((O,I,O), R),   4),
             case(Sphere((I,I,O), R),   8),
             case(Sphere((O,O,I), R),  16),
             case(Sphere((I,O,I), R),  32),
             case(Sphere((O,I,I), R),  64),
             case(Sphere((I,I,I), R), 128),
             // Should pick two voxels differing only in x-component
             case(CylinderX((O,O), R),   3),
             case(CylinderX((O,I), R),  48),
             case(CylinderX((I,O), R),  12),
             case(CylinderX((I,I), R), 192),
             // Should pick two voxels differing only in y-component
             case(CylinderY((O,O), R),   5),
             case(CylinderY((O,I), R),  80),
             case(CylinderY((I,O), R),  10),
             case(CylinderY((I,I), R), 160),
             // Should pick two voxels differing only in z-component
             case(CylinderZ((O,O), R),  17),
             case(CylinderZ((O,I), R),  68),
             case(CylinderZ((I,O), R),  34),
             case(CylinderZ((I,I), R), 136),

    )]
    fn which_voxels(roi: ROI, expected: usize) {
        let data = vec![1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0];
        let (l, n) = (4.0, 2);
        let vbox = VoxelBox::new((l,l,l), (n,n,n));
        let image = Image::new(vbox, data);
        let pattern = image.values_inside_roi(roi)
            .iter().sum::<Intensity>()
            as usize;
        println!("pattern {} {} expected", pattern, expected);
        assert!(pattern == expected);

    }
}

// TODO stop reinventing this wheel
fn mean(data: &[Intensity]) -> Option<Intensity> {
    data.iter().cloned().reduce(|a, b| a+b).map(|s| s / data.len() as Intensity)
}

fn mu_and_sigma(data: &[Intensity]) -> Option<(Intensity, Intensity)> {
    let mu = mean(data)?;
    let sigma = data.iter().cloned()
        .map(|x| x-mu)
        .map(|x| x*x)
        .sum::<Intensity>() / data.len() as Intensity;
    Some((mu, sigma))

}

#[cfg(test)]
mod test_mean {
    use super::*;

    #[test]
    fn test_mean() {
        let it = mean(&vec![1.0, 1.0, 1.0]);
        assert!(it == Some(1.0));

        let it = mean(&vec![1.0, 2.0, 3.0, 4.0]);
        assert!(it == Some(2.5));

        let it = mean(&vec![]);
        assert!(it == None);
    }

    #[test]
    fn test_mu_and_sigma() {
        if let Some((mu, sigma)) = mu_and_sigma(&vec![1.0, 1.0, 1.0]) {
            assert!(mu    == 1.0);
            assert!(sigma == 0.0);
        };

        if let Some((mu, sigma)) = mu_and_sigma(&vec![2.0, 4.0, 2.0, 4.0]) {
            assert!(mu    == 3.0);
            assert!(sigma == 1.0);
        };

        assert!(mu_and_sigma(&vec![]) == None);
    }
}

#[derive(Debug)]
pub struct FomConfig {
    pub rois: Vec<(ROI, Intensity)>,
    pub background_rois: Vec<ROI>,
    pub background_activity: Intensity
}

impl FomConfig {
    pub fn new(rois: Vec<(ROI, Intensity)>, background_rois: Vec<ROI>, background_activity: Intensity) -> Self {
        Self{ rois, background_rois, background_activity }
    }
}

#[allow(clippy::upper_case_acronyms)]
pub struct FOMS {
    pub crcs: Vec<Ratio>,
    pub snrs: Vec<Ratio>,
}

impl Image {
    pub fn foms(&self, config: &FomConfig, quiet: bool) -> FOMS {
        let FomConfig{ rois, background_rois, background_activity} = config;
        let background_measured = background_rois.iter().cloned()
            .map(|roi| mean(&self.values_inside_roi(roi)).unwrap())
            .sum::<Intensity>() / background_rois.len() as Intensity;
        if !quiet {println!("    background measured (set): {:4.1} ({})", background_measured, background_activity);}

        let mut crcs = vec![];
        let mut snrs = vec![];
        if !quiet {print!("    measured roi activities:");}
        for (roi, roi_activity) in rois.iter().cloned() {
            let (roi_measured, roi_sigma) = mu_and_sigma(&self.values_inside_roi(roi)).unwrap();
            if !quiet {print!("{:9.1} ({})", roi_measured, roi_activity);}
            crcs.push(100.0 *
                if roi_activity > *background_activity { // hot CRC
                    ((roi_measured / background_measured) - 1.0) /
                    ((roi_activity / background_activity) - 1.0)
                } else { // cold CRC
                    1.0 - (roi_measured / background_measured)
                }
            );

            snrs.push((roi_measured - background_measured) / roi_sigma); // TODO doesn't quite match antea
        }
        if !quiet {println!();}
        FOMS{crcs, snrs}
    }
}

#[cfg(test)]
mod test_crc {
    use super::*;

}
