use kiss3d;
use nalgebra as na;

use kiss3d::light::Light;
use kiss3d::window::Window;
use na::{Point3, Translation3};

use crate::weights as pet;

// TODO: x-axis is flipped on the display: fix it!

pub fn lor_weights(p1: pet::Point, p2: pet::Point, vbox: pet::VoxelBox) {
    let p13 = Point3::new(p1.x as f32, p1.y as f32, -1.0001);
    let p23 = Point3::new(p2.x as f32, p2.y as f32, -1.0001);

    let active_voxels = pet::WeightsAlongLOR::new(p1, p2, vbox.clone())
        .map(|(i, w)| ((i.x, i.y), w))
        .collect::<std::collections::HashMap<(usize, usize), pet::Length>>();

    let max_weight = active_voxels
        .iter()
        .map(|(_, w)| w)
        .max_by(|a,b| a.partial_cmp(b).expect("Weights contained NaN"))
        .unwrap_or(&1.0);

    let lor_colour = Point3::new(0.0, 1.0, 0.0);

    let mut window = Window::new("LOR weights");
    window.set_light(Light::StickToCamera);

    let mut voxels = Vec::with_capacity(vbox.n.x * vbox.n.y);
    {
        let vsize = vbox.voxel_size();
        let bsize = vbox.aabb.half_extents;
        let (bdx, bdy) = (bsize.x as f32, bsize.y as f32);
        let (vdx, vdy) = (vsize.x as f32, vsize.y as f32);
        let half_vbox = Translation3::new(-bdx, -bdy, 0.0);
        let half_voxel = Translation3::new(vdx / 2.0,
                                           vdy / 2.0,
                                           0.0);
        for x in 0..vbox.n.x {
            for y in 0..vbox.n.y {
                let mut v = window.add_cube(vdx * 0.99, vdy * 0.99, 1.0);
                v.append_translation(&half_vbox);
                v.append_translation(&half_voxel);
                v.append_translation(&Translation3::new(x as f32 * vdx, y as f32 * vdy, 0.0));
                if let Some(weight) = active_voxels.get(&(x,y)) {
                    v.set_color((*weight / max_weight) as f32, 0.0, 0.0);
                } else {
                    if x >= vbox.n.x || y >= vbox.n.y {
                    v.set_color(0.9, 0.9, 0.2);
                    } else {
                    v.set_color(0.2, 0.2, 0.2);
                    }
                }
                voxels.push(v);
            }
        }
    }

    let x_axis_lo = Point3::new(-90.0,   0.0, 0.0);
    let x_axis_hi = Point3::new( 90.0,   0.0, 0.0);
    let y_axis_lo = Point3::new(  0.0, -90.0, 0.0);
    let y_axis_hi = Point3::new(  0.0,  90.0, 0.0);
    let axis_colour = Point3::new(0.0, 0.0, 1.0);
    while window.render() {
        window.draw_line (&x_axis_lo, &x_axis_hi, &  axis_colour);
        window.draw_line (&y_axis_lo, &y_axis_hi, &  axis_colour);
        window.draw_line (&p13, &p23, &  lor_colour);
    }
}
