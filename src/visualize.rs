use kiss3d;
use nalgebra as na;

use kiss3d::light::Light;
use kiss3d::window::Window;
use na::{Point3, Translation3};

use crate::weights as pet;


pub fn lor_weights(p1: pet::Point, p2: pet::Point, vbox: pet::VoxelBox) {

    let active_voxels = pet::WeightsAlongLOR::new(p1, p2, vbox.clone())
        .map(|(i, w)| ((i.x, i.y, i.z), w))
        .collect::<std::collections::HashMap<(usize, usize, usize), pet::Length>>();

    let max_weight = active_voxels
        .iter()
        .map(|(_, w)| w)
        .max_by(|a,b| a.partial_cmp(b).expect("Weights contained NaN"))
        .unwrap_or(&1.0);

    let lor_colour = Point3::new(0.0, 1.0, 0.0);

    let mut window = Window::new("LOR weights");
    window.set_light(Light::StickToCamera);


    let vsize = vbox.voxel_size();
    let bsize = vbox.aabb.half_extents;
    let (bdx, bdy, bdz) = (bsize.x as f32, bsize.y as f32, bsize.z as f32);
    let (vdx, vdy, vdz) = (vsize.x as f32, vsize.y as f32, vsize.z as f32);
    let mut voxels = Vec::with_capacity(vbox.n.x * vbox.n.y);
    {
        let half_vbox = Translation3::new(-bdx, -bdy, -bdz);
        let half_voxel = Translation3::new(vdx / 2.0,
                                           vdy / 2.0,
                                           vdz / 2.0);
        let s = 0.99;
        for x in 0..vbox.n.x {
            for y in 0..vbox.n.y {
                for z in 0..vbox.n.z {
                    if let Some(weight) = active_voxels.get(&(x,y,z)) {
                        let mut v = window.add_cube(vdx * s, vdy * s, vdy * s);
                        v.append_translation(&half_vbox);
                        v.append_translation(&half_voxel);
                        v.append_translation(&Translation3::new(x as f32 * vdx, y as f32 * vdy, z as f32 * vdz));
                        v.set_color((*weight / max_weight) as f32, 0.0, 0.0);
                        voxels.push(v);
                        //v.set_material(material);
                    } else {
                        //v.set_color(0.1, 0.1, 0.1);
                    }
                }
            }
        }
    }

    let x_axis_lo = Point3::new(-bdx,   0.0, 0.0) * 1.5;
    let x_axis_hi = Point3::new( bdx,   0.0, 0.0) * 1.5;
    let y_axis_lo = Point3::new(  0.0, -bdy, 0.0) * 1.5;
    let y_axis_hi = Point3::new(  0.0,  bdy, 0.0) * 1.5;
    let axis_colour = Point3::new(0.0, 0.0, 1.0);

    let biggest = vbox.aabb.half_extents.x.max(vbox.aabb.half_extents.y) as f32;
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
