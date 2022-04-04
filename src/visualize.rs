use kiss3d::light::Light;
use kiss3d::window::Window;
use kiss3d::event::{Action, WindowEvent};
use kiss3d::scene::SceneNode;
use kiss3d::camera::ArcBall;
use kiss3d::nalgebra::{Point3, Translation3};

use crate::types::ns_to_ps;
use crate::types::Vectorf32;
use crate::types::{Time, Ratio};
use crate::weights::{FOV, LOR};

use structopt::clap::arg_enum;

arg_enum! {
    #[derive(Debug, Clone)]
    pub enum Shape {
        Box,
        Ball
    }
}

pub struct Scene {
    // Graphical state
    window: Window,
    voxels: Vec<SceneNode>,
    camera: ArcBall,
    // Helper
    //draw_lines: Box<dyn FnMut(&Window) -> ()>,
    // Parameters which define the scene
    lor: LOR,
    fov: FOV,
    lines: Vec<(Point3<f32>, Point3<f32>, Point3<f32>)>,
}

impl Scene {
    pub fn new(lor: LOR, fov: FOV) -> Self {
        let mut window = Window::new("LOR weights");
        window.set_light(Light::StickToCamera);

        // Axis lines
        let bsize = Vectorf32::from(fov.half_width);
        let (bdx, bdy, bdz) = (bsize.x as f32, bsize.y as f32, bsize.z as f32);
        let x_axis_lo = Point3::new(-bdx,   0.0, 0.0) * 1.5;
        let x_axis_hi = Point3::new( bdx,   0.0, 0.0) * 1.5;
        let y_axis_lo = Point3::new(  0.0, -bdy, 0.0) * 1.5;
        let y_axis_hi = Point3::new(  0.0,  bdy, 0.0) * 1.5;
        let z_axis_lo = Point3::new(  0.0,  0.0,-bdz) * 1.5;
        let z_axis_hi = Point3::new(  0.0,  0.0, bdz) * 1.5;
        let x_axis_colour = Point3::new(1.0, 0.0, 0.0);
        let y_axis_colour = Point3::new(0.0, 1.0, 0.0);
        let z_axis_colour = Point3::new(0.0, 0.0, 1.0);

        // LOR line
        let p1_f32 = Point3::new(lor.p1.x as f32, lor.p1.y as f32, lor.p1.z as f32);
        let p2_f32 = Point3::new(lor.p2.x as f32, lor.p2.y as f32, lor.p2.z as f32);
        let lor_colour = Point3::new(1.0, 1.0, 0.0);

        // FOV frame
        let w = Vectorf32::from(fov.half_width);
        let (bwx, bwy, bwz) = (w.x as f32, w.y as f32, w.z as f32);
        let box_000 = Point3::new(-bwx, -bwy, -bwz);
        let box_001 = Point3::new(-bwx, -bwy,  bwz);
        let box_010 = Point3::new(-bwx,  bwy, -bwz);
        let box_011 = Point3::new(-bwx,  bwy,  bwz);
        let box_100 = Point3::new( bwx, -bwy, -bwz);
        let box_101 = Point3::new( bwx, -bwy,  bwz);
        let box_110 = Point3::new( bwx,  bwy, -bwz);
        let box_111 = Point3::new( bwx,  bwy,  bwz);
        let box_colour = Point3::new(0.3, 0.3, 0.3);

        // Turn the above endpoints into actual lines
        let lines = vec![(x_axis_lo, x_axis_hi, x_axis_colour),
                         (y_axis_lo, y_axis_hi, y_axis_colour),
                         (z_axis_lo, z_axis_hi, z_axis_colour),
                         (p1_f32   , p2_f32  ,     lor_colour),
                         (box_000  , box_001 ,     box_colour),
                         (box_001  , box_011 ,     box_colour),
                         (box_011  , box_010 ,     box_colour),
                         (box_010  , box_000 ,     box_colour),
                         (box_100  , box_101 ,     box_colour),
                         (box_101  , box_111 ,     box_colour),
                         (box_111  , box_110 ,     box_colour),
                         (box_110  , box_100 ,     box_colour),
                         (box_000  , box_100 ,     box_colour),
                         (box_001  , box_101 ,     box_colour),
                         (box_011  , box_111 ,     box_colour),
                         (box_010  , box_110 ,     box_colour),

        ];

        Scene {
            window,
            voxels: vec![],
            camera: Self::init_camera(&fov),
            //draw_lines,
            lor,
            fov,
            lines,
        }
    }

    fn clear(&mut self) {
        println!("Clearing ... supposedly");
        for v in &mut self.voxels {
            self.window.remove_node(v);
            v.unlink();
        }
    }

    pub fn place_voxels(&mut self, shape: Shape, cutoff: Option<Ratio>, sigma: Option<Time>) {

        let active_voxels = self.lor.active_voxels(&self.fov, cutoff, sigma);

        let &max_weight = active_voxels
            .iter()
            .map(|(_, w)| w)
            .max_by(|a,b| a.partial_cmp(b).expect("Weights contained NaN"))
            .unwrap_or(&1.0);

        let vsize = Vectorf32::from(self.fov.voxel_size);
        let bsize = Vectorf32::from(self.fov.half_width);
        let (bdx, bdy, bdz) = (bsize.x as f32, bsize.y as f32, bsize.z as f32);
        let (vdx, vdy, vdz) = (vsize.x as f32, vsize.y as f32, vsize.z as f32);
        let mut voxels = Vec::with_capacity(self.fov.n[0] * self.fov.n[1]);
        let half_fov = Translation3::new(-bdx, -bdy, -bdz);
        let half_voxel = Translation3::new(vdx / 2.0,
                                           vdy / 2.0,
                                           vdz / 2.0);

        // Add voxel representations to the scene
        let s = 0.99;
        for (i, weight) in active_voxels {
            let relative_weight = (weight / max_weight) as f32;
            let mut v = match shape {
                Shape::Box  => self.window.add_cube(vdx * s, vdy * s, vdz * s),
                Shape::Ball => self.window.add_sphere(vdx.min(vdy).min(vdz) * relative_weight * 0.5),
            };
            v.append_translation(&half_fov);
            v.append_translation(&half_voxel);
            v.append_translation(&Translation3::new(i[0] as f32 * vdx, i[1] as f32 * vdy, i[2] as f32 * vdz));
            v.set_color(relative_weight, 0.1, 0.0);
            voxels.push(v);
            //v.set_material(material);
        }
    }

    fn draw_lines(&mut self) {
        for line in &self.lines {
            self.window.draw_line(&line.0, &line.1, &line.2);
        }
    }

    fn init_camera(fov: &FOV) -> ArcBall {
        // Fit image into FOV
        let half_width = Vectorf32::from(fov.half_width);
        let biggest = half_width.x.max(half_width.y) as f32;
        let distance = 4.0 * biggest;
        let zfar  = distance * 2.0;
        let znear = 0.1; //distance / 2.0;
        let fov = (biggest / distance).atan() * 2.2_f32;
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
                        println!("TODO: Toggle / change cutoff");
                    },
                    _ => {}
                }
            }
        }
    }
}

pub fn lor_weights(lor: LOR, fov: FOV, shape: Shape, cutoff: Option<Ratio>, sigma: Option<Time>) {
    let mut scene = Scene::new(lor, fov);
    scene.place_voxels(shape, cutoff, sigma);
    scene.main_loop();
}

pub fn vislor_command(fov: &FOV, lor: &LOR) -> String {
    let fov_half_width = Vectorf32::from(fov.half_width);
    format!(
        "cargo run --bin vislor -- box --fov-size {vx},{vy},{vz} --nvoxels {nx},{ny},{nz} --lor '{t1} {t2}   {x1} {y1} {z1}    {x2} {y2} {z2}'",
        vx = fov_half_width.x * 2.0,
        vy = fov_half_width.y * 2.0,
        vz = fov_half_width.z * 2.0,
        nx = fov.n[0],
        ny = fov.n[1],
        nz = fov.n[2],
        t1 = 0.0,
        t2 = ns_to_ps(lor.dt),
        x1 = lor.p1.x,
        y1 = lor.p1.y,
        z1 = lor.p1.z,
        x2 = lor.p2.x,
        y2 = lor.p2.y,
        z2 = lor.p2.z,
    )
}
