use structopt::StructOpt;
use structopt::clap::arg_enum;

use kiss3d;
use nalgebra as na;

use kiss3d::light::Light;
use kiss3d::window::Window;
use kiss3d::event::{Action, WindowEvent};
use kiss3d::scene::{SceneNode};
use kiss3d::camera::{ArcBall};
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

struct Scene {
    // Graphical state
    window: Window,
    voxels: Vec<SceneNode>,
    camera: ArcBall,
    // Helper
    //draw_lines: Box<dyn FnMut(&Window) -> ()>,
    // Parameters which define the scene
    p1: Point,
    p2: Point,
    t1: Time,
    t2: Time,
    vbox: VoxelBox,
    lines: Vec<(Point3<f32>, Point3<f32>, Point3<f32>)>,
}

impl Scene {
    fn new(t1: Time, t2: Time, p1: Point, p2: Point, vbox: VoxelBox) -> Self {
        let mut window = Window::new("LOR weights");
        window.set_light(Light::StickToCamera);

        // Define axis lines
        let bsize = vbox.half_width;
        let (bdx, bdy, bdz) = (bsize.x as f32, bsize.y as f32, bsize.z as f32);
        let x_axis_lo = Point3::new(-bdx,   0.0, 0.0) * 1.5;
        let x_axis_hi = Point3::new( bdx,   0.0, 0.0) * 1.5;
        let y_axis_lo = Point3::new(  0.0, -bdy, 0.0) * 1.5;
        let y_axis_hi = Point3::new(  0.0,  bdy, 0.0) * 1.5;
        let z_axis_lo = Point3::new(  0.0,  0.0,-bdz) * 1.5;
        let z_axis_hi = Point3::new(  0.0,  0.0, bdz) * 1.5;
        let x_axis_colour = Point3::new(1.0, 0.0, 0.0);
        let y_axis_colour = Point3::new(1.0, 1.0, 0.0);
        let z_axis_colour = Point3::new(0.0, 0.0, 1.0);

        // Define LOR
        let p1_f32 = Point3::new(p1.x as f32, p1.y as f32, p1.z as f32);
        let p2_f32 = Point3::new(p2.x as f32, p2.y as f32, p2.z as f32);
        let lor_colour = Point3::new(0.0, 1.0, 0.0);

        let lines = vec![(x_axis_lo, x_axis_hi, x_axis_colour),
                         (y_axis_lo, y_axis_hi, y_axis_colour),
                         (z_axis_lo, z_axis_hi, z_axis_colour),
                         (p1_f32   , p2_f32  ,   lor_colour),];

        Scene {
            window,
            voxels: vec![],
            camera: Self::init_camera(&vbox),
            //draw_lines,
            p1, p2, t1, t2,
            vbox,
            lines,
        }
    }

    fn clear(&mut self) {
        println!("Clearing ... supposedly");
        for mut v in &mut self.voxels {
            self.window.remove_node(&mut v);
            v.unlink();
        }
    }

    fn place_voxels(&mut self, args: &Cli) {
        // Did the user ask for TOF refinement
        let tof = match args.sigma {
            Some(sigma) => pet::tof_gaussian(self.p1, self.t1, self.p2, self.t2, &self.vbox, sigma),
            None        => Box::new(|x| x),
        };

        // Did the user ask to ignore low-weight voxels
        let threshold: Box<dyn FnMut(&(Index, Length)) -> bool> = match args.threshold {
            Some(thresh) => Box::new(move |(_, w)| w > &thresh),
            None         => Box::new(     |  _   |  true      ),
        };

        let active_voxels = pet::WeightsAlongLOR::new(self.p1, self.p2, &self.vbox)
            .map(tof)
            .filter(threshold)
            .collect::<std::collections::HashMap<Index, Length>>();

        let &max_weight = active_voxels
            .iter()
            .map(|(_, w)| w)
            .max_by(|a,b| a.partial_cmp(b).expect("Weights contained NaN"))
            .unwrap_or(&1.0);

        let vsize = self.vbox.voxel_size;
        let bsize = self.vbox.half_width;
        let (bdx, bdy, bdz) = (bsize.x as f32, bsize.y as f32, bsize.z as f32);
        let (vdx, vdy, vdz) = (vsize.x as f32, vsize.y as f32, vsize.z as f32);
        let mut voxels = Vec::with_capacity(self.vbox.n.x * self.vbox.n.y);
        let half_vbox = Translation3::new(-bdx, -bdy, -bdz);
        let half_voxel = Translation3::new(vdx / 2.0,
                                           vdy / 2.0,
                                           vdz / 2.0);

        // Add voxel representations to the scene
        let s = 0.99;
        for (i, weight) in active_voxels {
            let relative_weight = (weight / max_weight) as f32;
            let mut v = match args.shape {
                Shape::Box  => self.window.add_cube(vdx * s, vdy * s, vdy * s),
                Shape::Ball => self.window.add_sphere(vdx.min(vdy).min(vdy) * relative_weight * 0.8),
            };
            v.append_translation(&half_vbox);
            v.append_translation(&half_voxel);
            v.append_translation(&Translation3::new(i.x as f32 * vdx, i.y as f32 * vdy, i.z as f32 * vdz));
            v.set_color(relative_weight, 0.0, 0.0);
            voxels.push(v);
            //v.set_material(material);
        }
    }

    fn draw_lines(&mut self) {
        for line in &self.lines {
            self.window.draw_line(&line.0, &line.1, &line.2);
        }
    }

    fn init_camera(vbox: &VoxelBox) -> ArcBall {
        // Fit image into FOV
        let biggest = vbox.half_width.x.max(vbox.half_width.y) as f32;
        let distance = 4.0 * biggest;
        let zfar  = distance * 2.0;
        let znear = 0.1; //distance / 2.0;
        let fov = (biggest / distance).atan() * 2.2 as f32;
        kiss3d::camera::ArcBall::new_with_frustrum(
            fov, znear, zfar,
            Point3::new(0.0, 0.0, distance), // eye
            Point3::new(0.0, 0.0,      0.0), // at
        )
    }

    fn main_loop(&mut self) {
        while self.window.render_with_camera(&mut self.camera) {

            // Draw axes and LOR
            self.draw_lines();

            // Deal with events
            for event in self.window.events().iter() {
                use kiss3d::event::Key;
                match event.value {
                    WindowEvent::Key(Key::B, Action::Press, _) => {
                        println!("TODO: Toggle Box / Ball");
                        self.clear()
                    },
                    WindowEvent::Key(Key::S, Action::Press, _) => {
                        println!("TODO: Toggle TOF / change sigma");
                    },
                    WindowEvent::Key(Key::T, Action::Press, _) => {
                        println!("TODO: Toggle / change threshold");
                    },
                    _ => {}
                }
            }
        }

    }
}

pub fn lor_weights(t1: Time, t2: Time, p1: Point, p2: Point, vbox: VoxelBox) {

    let args = Cli::from_args();

    let mut scene = Scene::new(t1, t2, p1, p2, vbox);
    scene.place_voxels(&args);
    scene.main_loop();
}
