use structopt::{StructOpt, clap::arg_enum};
use ordered_float::NotNan;

arg_enum! {
    #[derive(Debug, Clone)]
    pub enum Phantom {
        Nema7,
        Jaszczak,
    }
}

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "nema7_foms", about = "Calculate NEMA7 Figures of Merit")]
pub struct Cli {
    /// Image file to analyse
    pub input_file: String,

    /// Which phantom is being analysed.
    #[structopt(possible_values = &Phantom::variants(), case_insensitive = true)]
    phantom: Phantom,
}

// --------------------------------------------------------------------------------
use std::error::Error;
use petalo::{fom::{InRoiFn, PointValue}, io::raw::Image3D};
use petalo::types::{Length, Intensity};
use petalo::mlem::Image;
use petalo::fom;
use petalo::fom::{Sphere, ROI, centres_of_slices_closest_to};

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::from_args();
    let image = Image3D::read_from_file(&args.input_file)?;
    let image = Image::from(&image);

    match args.phantom {
        Phantom::Nema7    =>    nema7_foms(&image),
        Phantom::Jaszczak => jaszczak_foms(&image),
    }
}

fn nema7_foms(image: &Image) -> Result<(), Box<dyn Error>> {

    let z_voxel_size = image.vbox.voxel_size[2];
    let z_half_width = image.vbox.half_width[2];

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

    // Background ROIs should be placed at z-slices closest to 0, ±1, ±2 cm.
    // These are the z-positions of the centres of the nearest slices
    let background_zs = centres_of_slices_closest_to(&[-20.0, -10.0, 0.0, 10.0, 20.0],
                                                     z_half_width, z_voxel_size);

    let mut background_roi_centres = vec![];
    for z in &background_zs {
        for (cx, cy) in &background_xys {
            background_roi_centres.push((*cx,*cy,*z));
        }
    };

    // The 6 hot spheres (position, diameter, activity)
    let r = 114.4 / 2.0; // displacement of centre of sphere from centre of body
    let spheres = vec![
        (sphere(r, 1, 10.0, 4.0)),
        (sphere(r, 2, 13.0, 4.0)),
        (sphere(r, 3, 17.0, 4.0)),
        (sphere(r, 4, 22.0, 4.0)),
        (sphere(r, 5, 28.0, 4.0)),
        (sphere(r, 0, 37.0, 4.0)),
    ];
    let foreground_z = fom::centre_of_slice_closest_to(z_half_width, z_voxel_size)(0.0);

    // The background count of the largest sphere (37mm) is also needed later
    // for the Accuracy of Corrections calculation. Create a place to store it.
    let mut bg_37 = None;

    let sphere_rois = spheres.iter()
        .map(|Sphere {x,y,r,..}| ROI::DiscZ((*x,*y,foreground_z), *r) )
        .collect::<Vec<_>>();

    // Annotate each voxel value with its 3D position
    let all_voxels = image.values_with_positions();

    // Discard voxels which do not lie inside any of the ROIs
    let relevant_voxels = discard_irrelevant_voxels(&sphere_rois, &background_roi_centres, &all_voxels);

    // Calculate the contrasts and background variabilities
    println!("Sphere diameter / mm    contrast / %   variability of background");
    for sphere in spheres {
        let (c,v) = contrast_and_variability(sphere, foreground_z,
                                             &background_xys, &background_zs, 1.0,
                                             &relevant_voxels).unwrap();
        if sphere.r == 37.0/2.0 { bg_37 = Some(c) }
        println!("{:20.0} {:10.1} {:15.1} ", sphere.r * 2.0, c, v);
    }

    // --- 7.4.2 ---------------------------------------------------------------------
    // ignore slices which lie within 30mm of ends.
    let hi_limit = 70.0         - 30.0;
    let lo_limit = 70.0 - 180.0 + 30.0;
    // Find voxel z-centres nearest to the limits
    let nearest = fom::centre_of_slice_closest_to(z_half_width, z_voxel_size);
    let hi_centre = nearest(hi_limit);
    let lo_centre = nearest(lo_limit);
    // If centre is beyond the limit, move inwards by one voxel
    let hi = if hi_centre <= hi_limit { hi_centre } else { nearest(hi_centre - z_voxel_size) };
    let lo = if lo_centre >= lo_limit { lo_centre } else { nearest(lo_centre + z_voxel_size) };
    // z-indices of first and last slices to be used
    let index_of = |z| fom::position_to_index(z, z_half_width, z_voxel_size);
    let   pos_of = |z| fom::index_to_position(z, z_half_width, z_voxel_size);
    let hi_index = index_of(hi);
    let lo_index = index_of(lo);

    // Ignore voxels which lie outside of the lung insert
    let lung = ROI::CylinderZ((0.0, 0.0), 30.0/2.0);
    let filter = lung.contains_fn();
    let lung_voxels: Vec<_> = fom::in_roi(filter, &all_voxels).collect();

    // For each z-slice divide mean within ROI, by 37mm background mean
    let bg_37 = bg_37.unwrap();
    let aocs = (lo_index..=hi_index).into_iter()
        .map(|i| { fom::mean_in_region(ROI::DiscZ((0.0, 0.0, pos_of(i)), 30.0/2.0), &lung_voxels) })
        .map(|v| 100.0 * v / bg_37)
        .collect::<Vec<_>>();

    println!("\nAOCs: {:.1?}", aocs);
    Ok(())
}


fn jaszczak_foms(image: &Image) -> Result<(), Box<dyn Error>> {

    let z_voxel_size = image.vbox.voxel_size[2];
    let z_half_width = image.vbox.half_width[2];

    // x-y centres of the background ROIs
    let background_xys = vec![
        ( 73.4, -32.0),
        ( 67.4,  35.2),
        (  1.5,  77.8),
        (-65.0,  38.2),
        (-69.6, -32.4),
        (- 2.4, -78.4),
        (  0.0,   0.0),
    ];

    // Ensure that the background ROIs are z-aligned with some slice
    let background_zs = centres_of_slices_closest_to(&[10.0, 20.0, 30.0, 40.0, 50.0, 60.0],
                                                     z_half_width, z_voxel_size);


    let mut background_roi_centres = vec![];
    for z in &background_zs {
        for (cx, cy) in &background_xys {
            background_roi_centres.push((*cx,*cy,*z));
        }
    };

    // The 6 hot spheres (position, diameter, activity)
    let r = 54.0; // displacement of centre of sphere from centre of body
    let spheres = vec![
        (sphere(r, 0,  9.5, 4.0)),
        (sphere(r, 1, 12.7, 4.0)),
        (sphere(r, 2, 15.9, 4.0)),
        (sphere(r, 3, 19.1, 4.0)),
        (sphere(r, 4, 25.4, 4.0)),
        (sphere(r, 5, 31.8, 4.0)),
    ];
    let foreground_z = fom::centre_of_slice_closest_to(z_half_width, z_voxel_size)(34.0);

    let sphere_rois = spheres.iter()
        .map(|Sphere {x,y,r,..}| ROI::DiscZ((*x,*y,foreground_z), *r) )
        .collect::<Vec<_>>();

    // Annotate each voxel value with its 3D position
    let all_voxels = image.values_with_positions();

    // Discard voxels which do not lie inside any of the ROIs
    let relevant_voxels = discard_irrelevant_voxels(&sphere_rois, &background_roi_centres, &all_voxels);

    // Calculate the contrasts and background variabilities
    println!("Sphere diameter / mm    contrast / %   variability of background / %");
    for sphere in spheres {
        let (c,v) = contrast_and_variability(sphere, foreground_z,
                                             &background_xys, &background_zs, 1.0,
                                             &relevant_voxels).unwrap();
        println!("{:20.1} {:10.1} {:15.1} ", sphere.r * 2.0, c, v);

    }
    Ok(())
}

/// Select voxels from the image which lie in any of the given ROIs, and
/// annotate them with their values
fn voxels_in_rois(image: Image, rois: Vec<ROI>) -> Vec<fom::PointValue> {
    // Annotate each voxel value with its 3D position
    let all_voxels = image.values_with_positions();
    //let mut relevant_voxels = vec![];
    todo!()
}

/// Place FOM sphere in Nth/6 angular position, with given diameter and activity
fn sphere(from_centre: Length, sphere_position: u16, diameter: Length, activity: Intensity) -> Sphere {
    let r = from_centre; // 114.4 / 2.0; // Radial displacement from centre
    let radians = std::f32::consts::TAU * sphere_position as f32;
    Sphere{x:r * radians.cos(), y:r * radians.sin(), r: diameter as Length / 2.0, a: activity}
}

fn contrast_and_variability(sphere: Sphere,
                            foreground_z: Length,
                            background_xys: &[(Length, Length)],
                            background_zs : &[Length],
                            bg_activity: f32,
                            voxels: &[PointValue],
) -> Option<(f32, f32)> {
    // Inspect single foreground ROI
    let Sphere { x, y, r, a: sphere_activity } = sphere;
    let sphere_mean = fom::mean_in_region(ROI::DiscZ((x, y, foreground_z), r), &voxels);
    // Inspect multiple background ROIs
    let mut bg_means = vec![];
    for (x,y) in background_xys {
        for z in background_zs {
            bg_means.push(fom::mean_in_region(ROI::DiscZ((*x, *y,*z), r), &voxels));
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

/// Make an f32 that is Ord. Panic if NaN
fn not_nan(f: f32) -> NotNan<f32> { NotNan::new(f).unwrap() }

/// Discard voxels which do not lie inside any of the ROIs
fn discard_irrelevant_voxels(
    sphere_rois           : &[ROI],
    background_roi_centres: &[(Length, Length, Length)],
    voxels                : &[PointValue],

) -> Vec<PointValue> {

    let max_roi_radius: Length = sphere_rois.iter()
        .map(ROI::r)
        .max_by_key(|&f| not_nan(f))
        .unwrap();

    // Background ROIs corresponding to biggest sphere
    let bg_rois = background_roi_centres.into_iter()
        .map(|&(x,y,z)| ROI::DiscZ((x,y,z), max_roi_radius))
        .collect::<Vec<_>>();

    // Individual functions checking whether point lies in a single background ROI
    let in_roi_checkers: Vec<InRoiFn> = bg_rois.iter()
        .chain(sphere_rois)
        .map(ROI::contains_fn)
        .collect::<Vec<_>>();

    // Combined function checking whether a point lies in any of the ROIs
    let in_some_roi = |p| in_roi_checkers.iter().any(|in_this_roi| in_this_roi(p));

    // Discard voxels which do not lie inside any of the ROIs
    voxels.iter().cloned()
        .filter(|(p, _)| in_some_roi(*p))
        .collect::<Vec<_>>()
}
