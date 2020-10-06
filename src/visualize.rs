use kiss3d;
use nalgebra as na;

use kiss3d::light::Light;
use kiss3d::window::Window;
use na::{Point3, Translation3};

use crate::weights as pet;
use crate::weights::{Length, Point, VoxelBox};


type N = f64;
fn make_gauss(mu: N, sigma: N) -> impl Fn(N) -> N {
    let root_two_pi = std::f64::consts::PI.sqrt();
    let a = 1.0 / (sigma * root_two_pi);
    move |x| {
        let y = (x - mu) / sigma;
        let z = y * y;
        a * (-0.5 * z).exp()
    }
}

type Time = f64;

type S = ((usize, usize, usize), Length);

const c : Length = 3e3; // cm / ns

fn tof(vbox: &pet::VoxelBox, t1: Time, t2: Time, p1: Point, p2: Point, mu: Length, sigma: Length) -> impl FnMut (&mut Option<Length>, S) -> Option<S> {
    let entry_point: Option<Point> = pet::entry(&p1, &p2, &vbox);
    let tof_centre_from_p1: Length = 0.5 * ((p1 - p2).norm() + c * (t1 - t2));
    let entry_point_from_p1: Option<Length> = entry_point.map(|ep| (ep-p1).norm());
    let original_distance_from_entry_to_tof_centre: Option<Length> = entry_point_from_p1.map(|ep| tof_centre_from_p1 - ep);
    let grr = original_distance_from_entry_to_tof_centre.unwrap();
    let gauss = make_gauss(mu, sigma);
    move |distance_from_tof_centre, (index, weight)| {
        let start: Length = match distance_from_tof_centre {
            None => grr,
            Some(d) => *d,
        };
        *distance_from_tof_centre = Some(start - weight);
        let midpoint: Length = start + weight / 2.0;
        Some((index, weight * gauss(midpoint)))
    }
}

pub fn lor_weights(t1: Time, t2: Time, p1: Point, p2: Point, vbox: VoxelBox) {

    let mu = 0.0;
    let sigma = 10.0;

    let active_voxels = pet::WeightsAlongLOR::new(p1, p2, &vbox)
        .map(|(i, w)| ((i.x, i.y, i.z), w))
        .scan(None, tof(&vbox, t1, t2, p1, p2, mu, sigma))
        .filter(|(_, w)| w > &0.01)
        .collect::<std::collections::HashMap<(usize, usize, usize), Length>>();

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
    let s = 0.99;
    for ((x,y,z), weight) in active_voxels {
        let relative_weight = (weight / max_weight) as f32;
        //let mut v = window.add_cube(vdx * s, vdy * s, vdy * s);
        let mut v = window.add_sphere(vdx.min(vdy).min(vdy) * relative_weight * 0.8);
        v.append_translation(&half_vbox);
        v.append_translation(&half_voxel);
        v.append_translation(&Translation3::new(x as f32 * vdx, y as f32 * vdy, z as f32 * vdz));
        v.set_color(relative_weight, 0.0, 0.0);
        voxels.push(v);
        //v.set_material(material);
    }

    let x_axis_lo = Point3::new(-bdx,   0.0, 0.0) * 1.5;
    let x_axis_hi = Point3::new( bdx,   0.0, 0.0) * 1.5;
    let y_axis_lo = Point3::new(  0.0, -bdy, 0.0) * 1.5;
    let y_axis_hi = Point3::new(  0.0,  bdy, 0.0) * 1.5;
    let axis_colour = Point3::new(0.0, 0.0, 1.0);

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

    let p1 = Point3::new(p1.x as f32, p1.y as f32, p1.z as f32);
    let p2 = Point3::new(p2.x as f32, p2.y as f32, p2.z as f32);

    while window.render_with_camera(&mut camera) {
        window.draw_line (&x_axis_lo, &x_axis_hi, &axis_colour);
        window.draw_line (&y_axis_lo, &y_axis_hi, &axis_colour);
        window.draw_line (&p1, &p2,               & lor_colour);
    }
}
