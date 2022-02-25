use std::path::PathBuf;
use structopt::StructOpt;
use petalo::utils::parse_range;
use petalo::io::hdf5::{Hdf5Lor, read_table};
use petalo::lorogram::{JustPhi, JustR, JustZ, Prompt, Scattergram, Lorogram};
use std::f32::consts::PI;


#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "show_logogram", about = "Interactive testing of logograms")]
pub struct Cli {

    #[structopt(short = "f", long)]
    pub input_file: PathBuf,

    /// The dataset location inside the input file
    #[structopt(short, long, default_value = "reco_info/lors")]
    pub dataset: String,

    /// Which rows of the input file should be loaded
    #[structopt(short, long, parse(try_from_str = parse_range::<usize>))]
    pub event_range: Option<std::ops::Range<usize>>,

}

fn main() -> Result<(), Box<dyn std::error::Error>> {

    let args = Cli::from_args();


    let infile  = args.input_file.into_os_string().into_string().unwrap();


    println!("===== z dependence ======================================");
    let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
    let sgram = fill_scattergram(JustZ::new(2000.0, 20), lors);

    for z in (-525..=525).step_by(50) {
        let p = (0.0, 0.0, z as f32);
        let v = sgram.value(p, p);
        println!("{z:5} {v:3.1}");
    }

    println!("===== phi dependence ====================================");
    let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
    let n = 15;
    let sgram = fill_scattergram(JustPhi::new(n), lors);

    for i in 0..n {
        let phi = PI * (i as f32 / n as f32);
        let x = phi.cos();
        let y = phi.sin();
        //println!("{:10.1} {:10.1} {:10.1} {:10.1} {}", x, y, phi, phi/PI, i);
        let p1 = (x, 0.0, 0.0);
        let p2 = (0.0, y, 0.0);
        let v = sgram.value(p1, p2);
        let phi = phi * 180.0 / PI;
        println!("{phi:5.1} {v:3.1}");
    }

    println!("===== r dependence ====================================");
    let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range)?;
    let nbins = 10;
    let r_max = 50.0;
    let sgram = fill_scattergram(JustR::new(r_max, nbins), lors);
    let step = r_max / nbins as f32;
    for i in 0..nbins {
        let r = (i as f32 + 0.5) * step;
        let p1 = (r,  100.0, 0.0);
        let p2 = (r, -100.0, 0.0);
        let v = sgram.value(p1, p2);
        println!("{r:5.1} {v:3.1}");
    }

    Ok(())
}

fn fill_scattergram<T: Lorogram + Clone>(prototype: T, lors: ndarray::Array1<Hdf5Lor>) ->  Scattergram<T> {
    let mut sgram = Scattergram::new(prototype);
    for Hdf5Lor { x1, y1, z1, x2, y2, z2, E1, E2, .. } in lors {
        if x1.is_nan() || x2.is_nan() { continue }
        let p1 = (x1, y1, z1);
        let p2 = (x2, y2, z2);
        let prompt = if E1.max(E2) < 511.0 { Prompt::Scatter } else { Prompt::True };
        sgram.fill(prompt, p1, p2);
    }
    sgram
}
