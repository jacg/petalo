use structopt::StructOpt;

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "nema7_analysis", about = "Calculate NEMA7 figures")]
pub struct Cli {
    /// Image file to analyse
    pub input_file: String,
}

// --------------------------------------------------------------------------------

use std::error::Error;
use petalo::io::raw::Image3D;
use petalo::types::{Length, Intensity};
use petalo::fom;

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let image = Image3D::read_from_file(&args.input_file)?;
    let image = petalo::mlem::Image::from(&image);

    let z_voxel_size = image.vbox.voxel_size[2];
    let z_half_width = image.vbox.half_width[2];

    // Background ROIs should be placed at z-slices closest to 0, ±1, ±2 cm.
    // These are the z-positions of the centres of the nearest slices
    let background_zs = centres_of_slices_closest_to(&[-20.0, -10.0, 0.0, 10.0, 20.0],
                                                     z_half_width, z_voxel_size);

    // Annotate each voxel value with its 3D position
    let pos_values = image.values_with_positions();

    // Collect voxels belonging to the 5 background ROI z-slices
    let mut slice = vec![vec![]; 5];
    for &(p,v) in &pos_values {
        for i in 0..5 {
            if p.z == background_zs[i] {
                slice[i].push((p,v));
                break;
            }
        }
    }

    // x-y centres of the 12 background ROIs
    let background_xys = vec![
        (   0.0, -86.0), // TODO all these should be 15 mm from edge of phantom
        (   0.0, -82.0),
        (  60.0, -82.0),
        ( -60.0, -82.0),
        (-100.0, -62.0),
        ( 100.0, -62.0),
        (-110.0, -20.0),
        ( 110.0, -20.0),
        ( -85.0,  40.0),
        (  85.0,  40.0),
        ( -52.0,  72.0),
        (  52.0,  72.0),
        (   0.0,  82.0),
    ];

    // The 6 hot spheres (position, diameter, activity)
    let spheres = vec![
        (nema7_sphere(1, 10, 4.0)),
        (nema7_sphere(2, 13, 4.0)),
        (nema7_sphere(3, 17, 4.0)),
        (nema7_sphere(4, 22, 4.0)),
        (nema7_sphere(5, 28, 4.0)),
        (nema7_sphere(0, 37, 4.0)),
    ];

    // The background count of the largest sphere (37mm) is also needed later
    // for the Accuracy of Corrections calculation. Create a place to store it.
    let mut bg_37 = None;

    // Calculate the contrasts and background variabilities
    println!("Sphere diameter / mm    contrast / %   variability of background");
    for sphere in spheres {
        let (c,v) = contrast_and_variability(sphere, &background_xys, &background_zs, &slice, 1.0).unwrap();
        if sphere.r == 37.0/2.0 { bg_37 = Some(c) }
        println!("{:20.0} {:10.1} {:15.1} ", sphere.r * 2.0, c, v);
    }

    // --- 7.4.2 ---------------------------------------------------------------------
    // ignore slices which lie within 30mm of ends.
    let hi_limit = 70.0         - 30.0;
    let lo_limit = 70.0 - 180.0 + 30.0;
    // Find voxel z-centres nearest to the limits
    let nearest = centre_of_slice_closest_to(z_half_width, z_voxel_size);
    let hi_centre = nearest(hi_limit);
    let lo_centre = nearest(lo_limit);
    // If centre is beyond the limit, move inwards by one voxel
    let hi = if hi_centre <= hi_limit { hi_centre } else { nearest(hi_centre - z_voxel_size) };
    let lo = if lo_centre >= lo_limit { lo_centre } else { nearest(lo_centre + z_voxel_size) };
    // z-indices of first and last slices to be used
    let index_of = |z| position_to_index(z, z_half_width, z_voxel_size);
    let   pos_of = |z| index_to_position(z, z_half_width, z_voxel_size);
    let hi_index = index_of(hi);
    let lo_index = index_of(lo);

    // Ignore voxels which lie outside of the lung insert
    let lung = fom::ROI::CylinderZ((0.0, 0.0), 30.0/2.0);
    let filter = lung.contains_fn();
    let lung_voxels: Vec<_> = in_roi(filter, &pos_values).collect();

    // For each z-slice divide mean within ROI, by 37mm background mean
    let bg_37 = bg_37.unwrap();
    let aocs = (lo_index..=hi_index).into_iter()
        .map(|i| { mean_in_region(fom::ROI::DiscZ((0.0, 0.0, pos_of(i)), 30.0/2.0), &lung_voxels) })
        .map(|v| 100.0 * v / bg_37)
        .collect::<Vec<_>>();

    println!("\nAOCs: {:.1?}", aocs);
    Ok(())
}

/// Place FOM sphere in Nth/6 angular position, with given diameter and activity
fn nema7_sphere(sphere_position: u16, diameter: u16, activity: Intensity) -> Sphere {
    let r = 114.4 / 2.0; // Radial displacement from centre
    let radians = std::f32::consts::TAU * sphere_position as f32;
    Sphere{x:r * radians.cos(), y:r * radians.sin(), r: diameter as Length / 2.0, a: activity}
}

/// x,y,r of NEMA7 sphere
#[derive(Clone, Copy)]
struct Sphere {
    x: Length,
    y: Length,
    r: Length,
    a: Intensity,
}

/// Mean of values associated with the voxels contained in the region
fn mean_in_region(roi: fom::ROI, voxels: &[fom::PointValue]) -> f32 {
    let filter = roi.contains_fn();
    let values: Vec<_> = in_roi(filter, voxels).map(|(_,v)| v).collect();
    fom::mean(&values).unwrap()
}

fn contrast_and_variability(sphere: Sphere,
                            background_xys: &[(Length, Length)],
                            background_zs : &[Length],
                            slices: &[Vec<fom::PointValue>],
                            bg_activity: f32,
) -> Option<(f32, f32)> {
    // Inspect single foreground ROI
    let Sphere { x, y, r, a: sphere_activity } = sphere;
    let sphere_mean = mean_in_region(fom::ROI::DiscZ((x, y, background_zs[2]), r), &slices[2]);
    // Inspect multiple background ROIs
    let mut bg_means = vec![];
    for (x,y) in background_xys {
        for (z, slice) in background_zs.iter().zip(slices) {
            bg_means.push(mean_in_region(fom::ROI::DiscZ((*x, *y,*z), r), &slice));
        }
    }
    // Calculate background variability
    let bg_mean = fom::mean(&bg_means)?;
    let bg_sd   = fom::sd  (&bg_means)?;
    let background_variability = 100.0 * bg_sd / bg_mean;
    // Calculate contrast
    let contrast = fom::crc(sphere_mean, sphere_activity, bg_mean, bg_activity);
    Some((contrast, background_variability))
}

/// Iterator which filters out voxels that lie outside given ROI
fn in_roi(in_roi: fom::InRoiFn, voxels: &[fom::PointValue]) -> impl Iterator<Item = fom::PointValue> + '_ {
    voxels.iter()
        .filter(move |(p,_)| in_roi(*p))
        .copied()
}

/// Convert 1D position to 1D index of containing voxel
fn position_to_index(position: Length, half_width: Length, voxel_size: Length) -> usize {
    ((position + half_width) / voxel_size) as usize
}

/// Convert 1D voxel index to 1D position of voxel's centre
fn index_to_position(index: usize, half_width: Length, voxel_size: Length) -> Length {
    (index as f32 + 0.5) * voxel_size - half_width
}

/// Return function which finds centre of nearest slice
fn centre_of_slice_closest_to(half_width: Length, voxel_size: Length) -> impl Fn(Length) -> Length {
    move |x| {
        let i = position_to_index(x, half_width, voxel_size);
        let x = index_to_position(i, half_width, voxel_size);
        x
    }
}

/// Adjust collection of 1D positions, to the centres of the nearest slices
fn centres_of_slices_closest_to(targets: &[Length], half_width: Length, voxel_size: Length) -> Vec<Length> {
    targets.iter().copied()
        .map(centre_of_slice_closest_to(half_width, voxel_size))
        .collect()
}
