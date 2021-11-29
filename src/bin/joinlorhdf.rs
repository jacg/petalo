use structopt::StructOpt;
use petalo::io;
use petalo::io::hdf5::Hdf5Lor;

#[derive(StructOpt, Debug, Clone)]
#[structopt(setting = structopt::clap::AppSettings::ColoredHelp)]
#[structopt(name = "joinlorhdf", about = "Combine datasets from separate HDF5 files into a single one")]
pub struct Cli {
    /// HDF5 input file with (t1 t2 x1 y1 z1 x2 y2 z2)
    pub inputs: Vec<String>,

    /// HDF5 output file for LORs found in input file
    #[structopt(short, long)]
    pub outfile: String,

    /// Group containing the dataset in each file
    #[structopt(short, long, default_value = "reco_info")]
    pub group: String,

    /// Dataset to read from each input file
    #[structopt(short, long, default_value = "lors")]
    pub dataset: String,

    // TODO allow using different group/dataset in output
}


type Data = Hdf5Lor; // TODO: add CLI switches for selecting type


fn main() -> hdf5::Result<()> {
    let args = Cli::from_args();
    let mut joined = Vec::<Data>::new();

    // ----- read data from separate files -------------------------------------------
    for filename in args.inputs {
        println!("Reading data from {}", filename);
        let path = format!("{}/{}", args.group, args.dataset);
        let mut data = io::hdf5::read_table::<Data>(&filename, &path, None)?;
        joined.extend_from_slice(data.as_slice_mut().unwrap());
        // TODO ndarray 0.14 -> 0.15: breaks our code in hdf5
        // joined.extend_from_slice(data.into_slice());
    }

    // --- write combined data to single file ----------------------------------------
    let outname = args.outfile;
    println!("Writing data to {}", outname);
    hdf5::File::create(outname)?
        .create_group(&args.group)?
        .new_dataset_builder()
        .with_data(&joined)
        .create(args.dataset.as_str())?;

    Ok(())
}
