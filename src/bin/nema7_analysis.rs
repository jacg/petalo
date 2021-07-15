use structopt::StructOpt;

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "nema7_analysis", about = "Calculate NEMA7 figures")]
pub struct Cli {
    /// Image file to analyse
    pub input_file: String,
}

// --------------------------------------------------------------------------------
// TODO change f32 to Length etc.
use std::error::Error;
use petalo::io::raw::Image3D;
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

    // The 6 hot spheres (position, diameter)
    let spheres = vec![
        (sphere(0, 37)),
        (sphere(1, 10)),
        (sphere(2, 13)),
        (sphere(3, 17)),
        (sphere(4, 22)),
        (sphere(5, 28)),
    ];

    // The background count of the largest sphere (37mm) is also needed later
    // for the Accuracy of Corrections calculation. Create a place to store it.
    let mut bg_37 = None;

    // Calculate the contrasts and background variabilities
    for sphere in spheres {
        let (c,v) = contrast_and_variability(sphere, &background_xys, &background_zs, &slice, 4.0, 1.0).unwrap();
        if sphere.r == 37.0/2.0 { bg_37 = Some(c) }
        println!("{:5.1}  {:5.1} ({:.0})", c, v, sphere.r * 2.0);
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
        .map(|i| { mean_in_region(0.0, 0.0, pos_of(i), 30.0/2.0, &lung_voxels) })
        .map(|v| 100.0 * v / bg_37)
        .collect::<Vec<_>>();

    println!("{:.1?}", aocs);
    Ok(())
}

fn sphere(sphere_position: u16, diameter: u16) -> HotSphere {
    let r = 114.4 / 2.0;
    let radians = std::f32::consts::TAU * sphere_position as f32;
    HotSphere{x:r * radians.cos(), y:r * radians.sin(), r: diameter as f32 / 2.0}
}

#[derive(Clone, Copy)]
struct HotSphere {
    x: f32,
    y: f32,
    r: f32,
}

fn mean_in_region(x: f32, y: f32, z: f32, r: f32, voxels: &[fom::PointValue]) -> f32 {
    let roi = fom::ROI::DiscZ((x,y,z), r);
    let filter = roi.contains_fn();
    let values: Vec<_> = in_roi(filter, voxels).map(|(_,v)| v).collect();
    fom::mean(&values).unwrap()
}

fn contrast_and_variability(sphere: HotSphere,
                            background_xys: &[(f32, f32)],
                            background_zs : &[f32],
                            slices: &[Vec<fom::PointValue>],
                            sphere_activity: f32,
                            bg_activity: f32,
) -> Option<(f32, f32)> {
    let HotSphere { x, y, r } = sphere;
    let sphere_mean = mean_in_region(x, y, background_zs[2], r, &slices[2]);

    let mut bg_means = vec![];
    for (x,y) in background_xys {
        for (z, slice) in background_zs.iter().zip(slices) {
            bg_means.push(mean_in_region(*x, *y,*z, r, &slice));
        }
    }
    let bg_mean = fom::mean(&bg_means)?;
    let bg_sd   = fom::sd  (&bg_means)?;
    let contrast = 100.0 * ((sphere_mean / bg_mean) - 1.0) / ((sphere_activity / bg_activity) - 1.0);
    let background_variability = 100.0 * bg_sd / bg_mean;
    Some((contrast, background_variability))
}

fn in_roi(in_roi: fom::InRoiFn, voxels: &[fom::PointValue]) -> impl Iterator<Item = fom::PointValue> + '_ {
    voxels.iter()
        .filter(move |(p,_)| in_roi(*p))
        .copied()
}

fn position_to_index(position: f32, half_width: f32, voxel_size: f32) -> usize {
    ((position + half_width) / voxel_size) as usize
}

fn index_to_position(index: usize, half_width: f32, voxel_size: f32) -> f32 {
    (index as f32 + 0.5) * voxel_size - half_width
}

fn centre_of_slice_closest_to(half_width: f32, voxel_size: f32) -> impl Fn(f32) -> f32 {
    move |x| {
        let i = position_to_index(x, half_width, voxel_size);
        let x = index_to_position(i, half_width, voxel_size);
        x
    }
}

fn centres_of_slices_closest_to(targets: &[f32], half_width: f32, voxel_size: f32) -> Vec<f32> {
    targets.iter().copied()
        .map(centre_of_slice_closest_to(half_width, voxel_size))
        .collect()
}
