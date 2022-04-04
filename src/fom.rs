use crate::io::raw;
use crate::types::{Lengthf32, Pointf32, Intensityf32, Ratiof32};
use crate::mlem::{Image, ImageData};
use crate::weights::FOV;

type BoxErr<T> = Result<T, Box<dyn std::error::Error>>;

pub fn load_image(filename: &std::path::Path, fov: FOV) -> BoxErr<Image> {
    let data = raw::read(filename)?.collect::<Result<_,_>>()?;
    Ok(Image::new(fov, data)) // TODO: Upgrade Image::new from panic to Result
}


#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum ROI {
    Sphere((Lengthf32, Lengthf32, Lengthf32), Lengthf32),
    CylinderX((Lengthf32, Lengthf32), Lengthf32),
    CylinderY((Lengthf32, Lengthf32), Lengthf32),
    CylinderZ((Lengthf32, Lengthf32), Lengthf32),
    DiscZ((Lengthf32, Lengthf32, Lengthf32), Lengthf32),
}

pub type InRoiFn = Box<dyn Fn(Pointf32) -> bool>;

impl ROI {

    pub fn contains_fn(&self) -> InRoiFn {
        match *self {
            ROI::Sphere((cx, cy, cz), radius) => Box::new(move |p: Pointf32| {
                let (x,y,z) = (p.x - cx, p.y - cy, p.z - cz);
                x*x + y*y + z*z < radius * radius
            }),

            ROI::CylinderX((cy, cz), radius) => Box::new(move |p: Pointf32| {
                let (y, z) = (p.y - cy, p.z - cz);
                y*y + z*z < radius*radius
            }),

            ROI::CylinderY((cx, cz), radius) => Box::new(move |p: Pointf32| {
                let (x, z) = (p.x - cx, p.z - cz);
                x*x + z*z < radius*radius
            }),

            ROI::CylinderZ((cx, cy), radius) => Box::new(move |p: Pointf32| {
                let (x, y) = (p.x - cx, p.y - cy);
                x*x + y*y < radius*radius
            }),

            ROI::DiscZ((cx, cy, z), radius) => Box::new(move |p: Pointf32| {
                let (x, y) = (p.x - cx, p.y - cy);
                z == p.z && x*x + y*y < radius*radius
            }),
        }
    }

    pub fn r(&self) -> Lengthf32 {
        match *self {
            ROI::Sphere   (_,r) => r,
            ROI::CylinderX(_,r) => r,
            ROI::CylinderY(_,r) => r,
            ROI::CylinderZ(_,r) => r,
            ROI::DiscZ    (_,r) => r,
        }
    }

}

/// A 3D point with an associated value. Used to represent voxels
pub type PointValue = (Pointf32, Intensityf32);

// TODO replace vec with iterator in output
impl Image {

    pub fn values_inside_roi(&self, roi: ROI) -> ImageData {
        let mut out = vec![];
        let roi_contains = roi.contains_fn();
        for (index, value) in self.data.iter().copied().enumerate() {
            let p = self.fov.voxel_centre1(index);
            if roi_contains(p) { out.push(value) }
        }
        out
    }

    pub fn values_with_positions(&self) -> Vec<PointValue> {
        self.data.iter().copied()
            .enumerate()
            .map(|(index, value)| (self.fov.voxel_centre1(index), value))
            .collect()
    }

}

/// Mean of values associated with the voxels contained in the region
pub fn mean_in_region(roi: ROI, voxels: &[PointValue]) -> f32 {
    let filter = roi.contains_fn();
    let values: Vec<_> = in_roi(filter, voxels).map(|(_,v)| v).collect();
    mean(&values).unwrap()
}

/// Iterator which filters out voxels that lie outside given ROI
pub fn in_roi(in_roi: InRoiFn, voxels: &[PointValue]) -> impl Iterator<Item = PointValue> + '_ {
    voxels.iter()
        .filter(move |(p,_)| in_roi(*p))
        .copied()
}

/// Convert 1D position to 1D index of containing voxel
pub fn position_to_index(position: Lengthf32, half_width: Lengthf32, voxel_size: Lengthf32) -> usize {
    ((position + half_width) / voxel_size) as usize
}

/// Convert 1D voxel index to 1D position of voxel's centre
pub fn index_to_position(index: usize, half_width: Lengthf32, voxel_size: Lengthf32) -> Lengthf32 {
    (index as f32 + 0.5) * voxel_size - half_width
}

/// Return function which finds centre of nearest slice
pub fn centre_of_slice_closest_to(half_width: Lengthf32, voxel_size: Lengthf32) -> impl Fn(Lengthf32) -> Lengthf32 {
    move |x| {
        let i = position_to_index(x, half_width, voxel_size);
                index_to_position(i, half_width, voxel_size)
    }
}

/// Adjust collection of 1D positions, to the centres of the nearest slices
pub fn centres_of_slices_closest_to(targets: &[Lengthf32], half_width: Lengthf32, voxel_size: Lengthf32) -> Vec<Lengthf32> {
    targets.iter().copied()
        .map(centre_of_slice_closest_to(half_width, voxel_size))
        .collect()
}

#[cfg(test)]
mod test_in_roi {
    use super::*;
    use rstest::rstest;

    // Arrange for outer voxels to be centred at +/- 100.0, when n = 10
    const MAGIC: Lengthf32 = 10.0 / 9.0 * 200.0;

    #[rstest(/**/    l ,  n,        centre        ,    r , expected_len,
             case(MAGIC, 10, (  0.0,   0.0,   0.0), 173.3, 1000), // r > sqrt(3) * 100; all voxel centres inside sphere
             case(MAGIC, 10, (  0.0,   0.0,   0.0), 173.2,  992), // r < sqrt(3) * 100; 992 = 8 corners missing
             case(MAGIC,  9, (  0.0,   0.0,   0.0), 173.2,  729), // all 9^3 included: coarser grid => outer centres closer
             case(MAGIC, 10, (200.0, 200.0, 200.0), 173.3,    1), // single voxel at one corner of box
             case(MAGIC, 10, (200.0, 200.0, 200.0), 173.2,    0), // slightly smaller r excludes the corner
    )]
    fn number_of_included_voxels_in_sphere(
        l: Lengthf32,
        n: usize,
        centre: (Lengthf32, Lengthf32, Lengthf32),
        r: Lengthf32,
        expected_len: usize
    ) {
        let data = vec![1.0; n*n*n];
        let fov = FOV::new((l,l,l), (n,n,n));
        let image = Image::new(fov, data);
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
        l: Lengthf32,
        n: usize,
        centre: (Lengthf32, Lengthf32),
        r: Lengthf32,
        expected_len: usize
    ) {
        let data = vec![1.0; n*n*n];
        let fov = FOV::new((l,l,l), (n,n,n));
        let image = Image::new(fov, data);
        let inside = image.values_inside_roi(ROI::CylinderZ(centre, r));
        println!("{} {}", inside.len(), expected_len);
        assert!(inside.len() == expected_len);
    }

    const I: Intensityf32 =  1.0;
    const O: Intensityf32 = -1.0;
    const R: Lengthf32 = 1.0;
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
        let fov = FOV::new((l,l,l), (n,n,n));
        let image = Image::new(fov, data);
        let pattern = image.values_inside_roi(roi)
            .iter().sum::<Intensityf32>()
            as usize;
        println!("pattern {} {} expected", pattern, expected);
        assert!(pattern == expected);

    }
}

// TODO stop reinventing this wheel
pub fn mean(data: &[Intensityf32]) -> Option<Intensityf32> {
    data.iter().cloned().reduce(|a, b| a+b).map(|s| s / data.len() as Intensityf32)
}

pub fn sd(data: &[Intensityf32]) -> Option<Intensityf32> {
    if data.len() < 2 { return None; }
    let mu = mean(data)?;
    let sum_of_deltas: Intensityf32 = data.iter()
        .map(|x| {let d = x-mu; d*d})
        .sum();
    Some(sum_of_deltas / data.len() as Intensityf32)
}

pub fn mu_and_sigma(data: &[Intensityf32]) -> Option<(Intensityf32, Intensityf32)> {
    let mu = mean(data)?;
    let sigma = data.iter().cloned()
        .map(|x| x-mu)
        .map(|x| x*x)
        .sum::<Intensityf32>() / data.len() as Intensityf32;
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

/// x,y,r of FOM sphere
#[derive(Clone, Copy)]
pub struct Sphere {
    pub x: Lengthf32,
    pub y: Lengthf32,
    pub r: Lengthf32,
    pub a: Intensityf32,
}

#[derive(Debug)]
pub struct FomConfig {
    pub rois: Vec<(ROI, Intensityf32)>,
    pub background_rois: Vec<ROI>,
    pub background_activity: Intensityf32
}

impl FomConfig {
    pub fn new(rois: Vec<(ROI, Intensityf32)>, background_rois: Vec<ROI>, background_activity: Intensityf32) -> Self {
        Self{ rois, background_rois, background_activity }
    }
}

#[allow(clippy::upper_case_acronyms)]
pub struct FOMS {
    pub crcs: Vec<Ratiof32>,
    pub snrs: Vec<Ratiof32>,
}

#[allow(clippy::upper_case_acronyms)]
pub struct FOM {
    pub r: Lengthf32,
    pub crc: Ratiof32,
    pub bg_variability: Ratiof32,
    pub snr: Ratiof32,
}

impl Image {

    pub fn foms(&self, config: &FomConfig, quiet: bool) -> FOMS {
        let FomConfig{ rois, background_rois, background_activity} = config;
        let background_measured = background_rois.iter().cloned()
            .map(|roi| mean(&self.values_inside_roi(roi)).unwrap())
            .sum::<Intensityf32>() / background_rois.len() as Intensityf32;
        if !quiet {println!("    background measured (set): {:4.1} ({})", background_measured, background_activity);}

        let mut crcs = vec![];
        let mut snrs = vec![];
        if !quiet {print!("    measured roi activities:");}
        for (roi, roi_activity) in rois.iter().cloned() {
            let (roi_measured, roi_sigma) = mu_and_sigma(&self.values_inside_roi(roi)).unwrap();
            if !quiet {print!("{:9.1} ({})", roi_measured, roi_activity);}
            crcs.push(crc(roi_measured, roi_activity, background_measured, *background_activity));
            snrs.push((roi_measured - background_measured) / roi_sigma); // doesn't quite match antea
        }
        if !quiet {println!();}
        FOMS{crcs, snrs}
    }

}

/// Calculate hot or cold Contrast Recovery Coefficient as percentage.
///
/// The exact calculation performed depends on whether the ROI is formally
/// hotter or colder than the background.
pub fn crc(roi_measured: Intensityf32, roi_activity: Intensityf32,
           bgd_measured: Intensityf32, bgd_activity: Intensityf32) -> Ratiof32 {
    100.0 * if roi_activity > bgd_activity { // hot CRC
        ((roi_measured / bgd_measured) - 1.0) /
        ((roi_activity / bgd_activity) - 1.0)
    } else { // cold CRC
        1.0 - (roi_measured / bgd_measured)
    }
}

#[cfg(test)]
mod test_crc {
    //use super::*;

}
