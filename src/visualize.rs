use structopt::StructOpt;
use structopt::clap::arg_enum;

use kiss3d;
use nalgebra as na;

use kiss3d::light::Light;
use kiss3d::window::Window;
use na::{Point3, Translation3};

use crate::weights as pet;
use crate::weights::{Length, Time, Point, VoxelBox, Index};

arg_enum! {
    #[derive(Debug)]
    enum Shape {
        Box,
        Ball
    }
}

#[derive(StructOpt, Debug)]
#[structopt(name = "petalo", about = "3d PET reconstruction with TOF")]
pub struct Cli {

    /// Apply TOF correction with given sensitivity. If not sepcified, no TOF
    /// correction is performed.
    #[structopt(short, long)]
    sigma: Option<Length>,

    /// Ignore voxels with weight below this threshold.
    #[structopt(short, long)]
    threshold: Option<Length>,

    /// How to represent voxels. BOX is better for viewing the geometric
    /// weights; BALL is better for viewing TOF weights.
    #[structopt(possible_values = &Shape::variants(), case_insensitive = true, default_value = "box")]
    shape: Shape,

}


pub fn lor_weights(t1: Time, t2: Time, p1: Point, p2: Point, vbox: VoxelBox) {

    let args = Cli::from_args();

    // Did the user ask for TOF refinement
    let tof = match args.sigma {
        Some(sigma) => pet::tof_gaussian(p1, t1, p2, t2, &vbox, sigma),
        None        => Box::new(|x| x),
    };

    // Did the user ask to ignore low-weight voxels
    let threshold: Box<dyn FnMut(&(Index, Length)) -> bool> = match args.threshold {
        Some(thresh) => Box::new(move |(_, w)| w > &thresh),
        None         => Box::new(     |  _   |  true      ),
    };

    let active_voxels = pet::WeightsAlongLOR::new(p1, p2, &vbox)
        .map(tof)
        .filter(threshold)
        .collect::<std::collections::HashMap<Index, Length>>();

    let &max_weight = active_voxels
        .iter()
        .map(|(_, w)| w)
        .max_by(|a,b| a.partial_cmp(b).expect("Weights contained NaN"))
        .unwrap_or(&1.0);

    let lor_colour = Point3::new(0.0, 1.0, 0.0);

    let mut window = Window::new("LOR weights");
    window.set_light(Light::StickToCamera);


    let vsize = vbox.voxel_size;
    let bsize = vbox.half_width;
    let (bdx, bdy, bdz) = (bsize.x as f32, bsize.y as f32, bsize.z as f32);
    let (vdx, vdy, vdz) = (vsize.x as f32, vsize.y as f32, vsize.z as f32);
    let mut voxels = Vec::with_capacity(vbox.n.x * vbox.n.y);
    let half_vbox = Translation3::new(-bdx, -bdy, -bdz);
    let half_voxel = Translation3::new(vdx / 2.0,
                                       vdy / 2.0,
                                       vdz / 2.0);

    // Add voxel representations to the scene
    let s = 0.99;
    for (i, weight) in active_voxels {
        let relative_weight = (weight / max_weight) as f32;
        let mut v = match args.shape {
            Shape::Box  => window.add_cube(vdx * s, vdy * s, vdy * s),
            Shape::Ball => window.add_sphere(vdx.min(vdy).min(vdy) * relative_weight * 0.8),
        };
        v.append_translation(&half_vbox);
        v.append_translation(&half_voxel);
        v.append_translation(&Translation3::new(i.x as f32 * vdx, i.y as f32 * vdy, i.z as f32 * vdz));
        v.set_color(relative_weight, 0.0, 0.0);
        voxels.push(v);
        //v.set_material(material);
    }

    // Define axis lines
    let x_axis_lo = Point3::new(-bdx,   0.0, 0.0) * 1.5;
    let x_axis_hi = Point3::new( bdx,   0.0, 0.0) * 1.5;
    let y_axis_lo = Point3::new(  0.0, -bdy, 0.0) * 1.5;
    let y_axis_hi = Point3::new(  0.0,  bdy, 0.0) * 1.5;
    let axis_colour = Point3::new(0.0, 0.0, 1.0);

    // Fit image into FOV
    let biggest = vbox.half_width.x.max(vbox.half_width.y) as f32;
    let distance = 4.0 * biggest;
    let zfar  = distance * 2.0;
    let znear = 0.1; //distance / 2.0;
    let fov = (biggest / distance).atan() * 2.2 as f32;
    let mut camera = kiss3d::camera::ArcBall::new_with_frustrum(
        fov, znear, zfar,
        Point3::new(0.0, 0.0, distance), // eye
        Point3::new(0.0, 0.0,      0.0), // at
    );

    // Define LOR
    let p1 = Point3::new(p1.x as f32, p1.y as f32, p1.z as f32);
    let p2 = Point3::new(p2.x as f32, p2.y as f32, p2.z as f32);

    while window.render_with_camera(&mut camera) {

        // Draw axes and LOR
        window.draw_line (&x_axis_lo, &x_axis_hi, &axis_colour);
        window.draw_line (&y_axis_lo, &y_axis_hi, &axis_colour);
        window.draw_line (&p1, &p2,               & lor_colour);
    }
}
