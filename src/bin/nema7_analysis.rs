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
    let image_io   = Image3D::read_from_file(&args.input_file)?;
    let image_mlem = petalo::mlem::Image::from(&image_io);

    println!("{:?} {:?}", image_io.pixels, image_io.mm);
    println!("{:?}", image_mlem.vbox);

    let voxel_size = image_mlem.vbox.voxel_size[2];
    let half_width = image_mlem.vbox.half_width[2];
    let z_targets = vec![-20.0, -10.0, 0.0, 10.0, 20.0];
    let z_indices = z_targets.into_iter()
        .map(|z| position_to_index(z, half_width, voxel_size))
        .collect::<Vec<_>>();
    println!("{:?}", z_indices);
    let background_zs = z_indices.iter().copied()
        .map(|i| index_to_position(i, half_width, voxel_size))
        .collect::<Vec<_>>();
    println!("{:?}", background_zs);

    let z = background_zs[3];

    let pos_values = image_mlem.values_with_positions();
    let mut layer = vec![vec![]; 5];
    for (p,v) in pos_values {
        for i in 0..5 {
            if p.z == background_zs[i] {
                layer[i].push((p,v));
                break;
            }
        }
    }


    let background_xys = vec![
        (   0.0, -86.0),
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

    let r = 114.4 / 2.0;
    let xyr = |sphere_position, diameter| {
        let radians = std::f32::consts::TAU * sphere_position as f32;
        HotSphere{x:r * radians.cos(), y:r * radians.sin(), r: diameter as f32 / 2.0}
    };

    let hot_specs = vec![
        (xyr(0, 37)),
        (xyr(1, 10)),
        (xyr(2, 13)),
        (xyr(3, 17)),
        (xyr(4, 22)),
        (xyr(5, 28)),

    ];

    for sphere in hot_specs {
        let (c,v) = crc(sphere, &background_xys, &background_zs, &layer, 4.0, 1.0).unwrap();
        println!("{:5.1}  {:5.1} ({:.0})", c, v, sphere.r * 2.0);
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct HotSphere {
    x: f32,
    y: f32,
    r: f32,
}

fn crc(sphere: HotSphere,
       background_xys: &[(f32, f32)],
       background_zs : &[f32],
       layers: &[Vec<fom::PointValue>],
       sphere_activity: f32,
       bg_activity: f32,
) -> Option<(f32, f32)> {
    let HotSphere { x, y, r } = sphere;
    let sphere_roi = fom::ROI::DiscZ((x,y, background_zs[2]), r);
    let sphere_filter = sphere_roi.contains_fn();
    let sphere_values: Vec<_> = in_roi(sphere_filter, &layers[2]).map(|(_,v)| v).collect();
    let sphere_mean = fom::mean(&sphere_values)?;

    let mut means = vec![];

    for (x,y) in background_xys {
        for (z, layer) in background_zs.iter().zip(layers) {
            let bg_roi = fom::ROI::DiscZ((*x,*y,*z), r);
            let bg_filter = bg_roi.contains_fn();
            let bg_values: Vec<_> = in_roi(bg_filter, &layer).map(|(_,v)| v).collect();
            let bg_mean = fom::mean(&bg_values)?;
            means.push(bg_mean);
        }
    }
    let bg_mean = fom::mean(&means)?;
    let bg_sd   = fom::sd  (&means)?;
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
