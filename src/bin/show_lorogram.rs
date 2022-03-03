use std::path::PathBuf;
use structopt::StructOpt;
use petalo::utils::parse_range;
use petalo::weights::LOR;
use petalo::io::hdf5::{Hdf5Lor, read_table};
use petalo::lorogram::{JustDeltaZ, JustPhi, JustR, JustZ, ZAndDeltaZ, Prompt, Scattergram, Lorogram};
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

    {
        println!("===== z dependence ======================================");
        let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
        let nbins = 20;
        let l = 200.0;
        let sgram = fill_scattergram(&|| Box::new(JustZ::new(l, nbins)), lors);
        let l0 = - l / 2.0;
        let dl = l / (nbins as f32);

        println!("     z       (s/t) + 1     trues   scatters");
        for i in 0..nbins {
            let z = l0 + (i as f32 + 0.5) * dl;
            let p = (0.0, 0.0, z as f32);
            let (v, t, s) = sgram.triplet(&LOR::from((p, p)));
            println!("{z:7.1}   {v:10.2}    {t:8}  {s:8}");
        }
    }
    {
        println!("===== phi dependence ====================================");
        let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
        let nbins = 15;
        let sgram = fill_scattergram(&|| Box::new(JustPhi::new(nbins)), lors);

        println!("   phi       (s/t) + 1     trues   scatters");
        for i in 0..nbins {
            let phi = PI * ((i as f32 + 0.5) / nbins as f32);
            let x = phi.cos();
            let y = phi.sin();
            let p1 = (x, 0.0, 0.0);
            let p2 = (0.0, y, 0.0);
            let (v, t, s) = sgram.triplet(&LOR::from((p1, p2)));
            let phi_in_degrees = phi * 180.0 / PI;
            println!("{phi_in_degrees:7.1}   {v:10.2}    {t:8}  {s:8}");
        }
    }
    {
        println!("===== r dependence ====================================");
        let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
        let nbins = 15;
        let r_max = 120.0;
        let sgram = fill_scattergram(&|| Box::new(JustR::new(r_max, nbins)), lors);
        let step = r_max / nbins as f32;
        println!("     r       (s/t) + 1     trues   scatters");
        for i in 0..nbins {
            let r = (i as f32 + 0.5) * step;
            let p1 = (r,  100.0, 0.0);
            let p2 = (r, -100.0, 0.0);
            let (v, t, s) = sgram.triplet(&LOR::from((p1, p2)));
            println!("{r:7.1}   {v:10.2}    {t:8}  {s:8}");
        }
    }
    {
        println!("===== obliqueness ====================================");
        let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
        let nbins = 20;
        let dz_max = 1000.0;
        let sgram = fill_scattergram(&|| Box::new(JustDeltaZ::new(dz_max, nbins)), lors);
        let step = dz_max / nbins as f32;
        println!("     dz      (s/t) + 1     trues   scatters");
        for i in 0..nbins {
            let dz = (i as f32 + 0.5) * step;
            let p1 = (0.0, 0.0,  dz/2.0);
            let p2 = (0.0, 0.0, -dz/2.0);
            let (v, t, s) = sgram.triplet(&LOR::from((p1, p2)));
            println!("{dz:7.1}   {v:10.2}    {t:8}  {s:8}");
        }
    }
    {
        println!("===== z and dz ====================================");
        let lors = read_table::<Hdf5Lor>(&infile, &args.dataset, args.event_range.clone())?;
        let (nbins_z, nbins_dz) = (20, 20);
        let (l, dz_max) = (200.0, 1000.0);
        let sgram = fill_scattergram(&|| Box::new(ZAndDeltaZ::new(l, nbins_z, dz_max, nbins_dz)), lors);
        let l0 = - l / 2.0;
        let dl = l / (nbins_z as f32);
        let step = dz_max / nbins_dz as f32;
        for i in 0..nbins_z {
            let z = (l0 + (i as f32 + 0.5) * dl) as f32;
            print!("z={z:6.1}    ");
            for j in 0..nbins_dz {
                let dz = (j as f32 + 0.5) * step;
                let p1 = (0.0, 0.0, z + dz/2.0);
                let p2 = (0.0, 0.0, z - dz/2.0);
                let v = sgram.value(&LOR::from((p1, p2)));
                print!(" {v:6.1}");
            }
            println!();
        }
    }

    Ok(())
}

fn fill_scattergram(make_empty_lorogram: &(dyn Fn() -> Box<dyn Lorogram>), lors: ndarray::Array1<Hdf5Lor>) ->  Scattergram {
    let mut sgram = Scattergram::new(make_empty_lorogram);
    for Hdf5Lor { x1, y1, z1, x2, y2, z2, E1, E2, .. } in lors {
        if x1.is_nan() || x2.is_nan() { continue }
        let p1 = (x1, y1, z1);
        let p2 = (x2, y2, z2);
        let prompt = if E1.min(E2) < 511.0 { Prompt::Scatter } else { Prompt::True };
        sgram.fill(prompt, &LOR::from((p1, p2)));
    }
    sgram
}
