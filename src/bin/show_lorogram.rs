use std::path::PathBuf;
use structopt::StructOpt;
use petalo::utils::parse_range;
use petalo::io::hdf5::{Hdf5Lor, read_table};
use petalo::lorogram::{JustZ, Prompt, Scattergram};


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

    let prototype = JustZ::new(2000.0, 20);

    let infile  = args.input_file.into_os_string().into_string().unwrap();

    let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range)?;
    let mut sgram = Scattergram::new(prototype.clone());
    for Hdf5Lor { x1, y1, z1, x2, y2, z2, E1, E2, .. } in lors {
        if x1.is_nan() || x2.is_nan() { continue }
        let p1 = (x1, y1, z1);
        let p2 = (x2, y2, z2);
        let prompt = if E1.max(E2) < 511.0 { Prompt::Scatter } else { Prompt::True };
        sgram.fill(prompt, p1, p2);
    }

    for z in (-525..=525).step_by(50) {
        let p = (0.0, 0.0, z as f32);
        let v = sgram.value(p, p);
        println!("{v:3.1}");
    }

    Ok(())
}
