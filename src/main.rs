use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::path::PathBuf;

// ----- CLI --------------------------------------------------
use clap::Parser;

/// Try to parse SimSET (custom) history files
#[derive(Parser, Debug)]
struct Args {

    /// History file to parse
    file: PathBuf,

    /// Maximum number of decays to parse
    #[arg(short = 'n', long, default_value_t = 10)]
    max_decays: usize,

    #[arg(value_enum, short = 't', long = "file-type", default_value_t = HistoryFileType::Custom)]
    history_file_type: HistoryFileType,

}

#[derive(clap::ValueEnum, Clone, Debug)]
enum HistoryFileType {
    Standard,
    Custom,
}

impl Display for HistoryFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HistoryFileType::Standard => f.write_str("Standard"),
            HistoryFileType::Custom   => f.write_str("Custom"),
        }
    }
}

// ----- Parser --------------------------------------------------
use binrw::binrw;
use binrw::io::{Seek, SeekFrom};
use binrw::BinReaderExt;

#[binrw]
#[derive(Debug)]
enum Standard {

    #[br(magic = 1_u8)]
    Decay {
    pos: Point,
    weight: f64,
    time: f64,
    decay_type: u64, // Needs to be mapped to an enum (PhgEn_DecayTypeTy)
    },

    #[br(magic = 2_u8)]
    Photon { // Both blue and pink?
        bytes: [u8; 72],
    },
}

#[binrw]
#[derive(Debug)]
enum Custom {

    #[br(magic = 0_u8)]
    Uhh,

    #[br(magic = 1_u8)]
    Photon {
        pos: Point,
        // x-cos
        // y-cos
        // z-cos
        // scatters in object
        // scatters in collimator
        // decay_weight: f64,
        // weight: f64,
        energy: f64,
        // travel_distance: f64,
        // decay_pos: Point,
        decay_time: f64,
        // decay_type: u32,

        // LOR
        // transaxial_distance: f64
        // azimuthal angle index
        // axial_position: f64,

        // detector_pos: Point,
        // detector_angle: f64,
        // detector crystal
        num_detector_interactions: i32,
        num_detector_interactions_again: i32,
        #[br(count = num_detector_interactions)]
        detector_interactions: Vec<DetectorInteraction>,
    },

}

#[binrw]
#[derive(Debug)]
struct Point {
    x: f64,
    y: f64,
    z: f64,
}

#[binrw]
#[derive(Debug)]
struct DetectorInteraction {
    pos: Point,
    energy_deposited: f64,
    is_active: u8,
    zero: u8, // There seems to be a trailing zero
}

fn custom(args: Args) -> Result<(), Box<dyn Error>> {
    let mut x = File::open(args.file)?;
    x.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    while let Ok(y) = x.read_le::<Custom>() {
        use Custom::*;
        match y {
            Uhh => {
                if count >= args.max_decays { break }
                count += 1;
                println!("============================================================");
            },
            Photon { pos: Point { x, y, z }, energy, decay_time, num_detector_interactions: n, detector_interactions, .. } => {
                println!("\n   ({x:7.2} {y:7.2} {z:7.2})  E: {energy:6.2}    decay time: {decay_time:5.2}");
                for (n, DetectorInteraction { pos: Point { x, y, z }, energy_deposited: energy, .. }) in detector_interactions.iter().enumerate() {
                    println!("{:2} ({x:7.2} {y:7.2} {z:7.2})     {energy:6.2}", n+1);
                }
            },
        }
    }
    Ok(())

}

fn standard(args: Args) -> Result<(), Box<dyn Error>> {
    let mut x = File::open(args.file)?;
    x.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    while let Ok(y) = x.read_le::<Standard>() {
        use Standard::*;
        match y {
            Decay { pos: Point { x, y, z }, weight, time, decay_type } => {
                if count >= args.max_decays { break }
                count += 1;
                println!("--------------------------------------");
                println!("{x:5.2} {y:5.2} {z:5.2}   {weight:5.2}    {time:5.2}   {decay_type}");
            },
            Photon { bytes } => {
                println!("");
                for byte in bytes { print!("{byte:02x} "); }
                println!("");
                for offset in 0..4 {
                    let ints = offset_by_as(&bytes, offset);
                    show_ints(&ints, offset);
                }
                println!("");
                for offset in 0..8 {
                    let reals = offset_by_as(&bytes, offset);
                    show_reals(&reals, offset);
                }
                println!("");
                for offset in 0..4 {
                    let reals = offset_by_as(&bytes, offset);
                    show_reals_32(&reals, offset);
                }
            },
        }
    }
    Ok(())
}


// typedef struct  {
//     double x_position;
//     double y_position;
//     double z_position;
// } PHG_Position;

// typedef struct {
//     PHG_Position        location;       /* Origination of decay */
//     double              startWeight;    /* starting weight for this decay */
//     /* double           eventWeight;    old decay weight for possible changes during tracking */
//     double              decayTime;      /* time, in seconds, between scan start and the decay */
//     PhgEn_DecayTypeTy   decayType;      /* single photon/positron/multi-emission/random etc. */
// } PHG_Decay;

// typedef enum {
//     PhgEn_SinglePhoton,
//     PhgEn_Positron,
//     PhgEn_PETRandom,    /* Artificial random coincidence events */
//     PhgEn_Complex,      /* For isotopes with multiple possible decay products - not implemented */
//     PhgEn_Unknown       /* for unassigned or error situations */
// } PhgEn_DecayTypeTy;


// Skip the header and collect a kilobyte into the xxx file
//
// cat Results/Test/SimSET_Sim_Discovery_ST/division_0/phg_hf.hist | dd bs=1 skip=32768 count=1024 > xxx

// at the beginning of each record in the history file is a flag-byte which
// indicates if the following record is a decay, a blue photon, or a pink
// photon. Processing the history file requires reading a flag byte to determine
// what type of record follows and then reading the appropriate amount of data.
//
// cat xxx | od --endian=little --skip-bytes 0 --read-bytes 1 -tx1     ->  01
// 1 presumably means it's a decay.

// ... a PHG_decay starts with a PHG_position which is three doubles
// cat xxx | od --endian=little --skip-bytes 1 --read-bytes 24 -tfD    ->   0.03716783271826184 0.12195498170336516 -0.0670295889395485

// ... then there's double for the starting weight of the decay
// cat xxx | od --endian=little --skip-bytes 25 --read-bytes 8 -tfD    ->   1.0000000016947894

// ... and another double for the decay time
// cat xxx | od --endian=little --skip-bytes 33 --read-bytes 8 -tfD    ->   18.77595793825161

// Next  a PhgEn_DecayTypeTy enum, so let's guess it's a byte (after further inspection, it looks more like a little endian 8-byte int)
// cat xxx | od --endian=little --skip-bytes 41 --read-bytes 1 -tx1    ->   01
// which looks like it's a positron

//

fn offset_by_as<T: Clone>(bytes: &[u8], offset: usize) -> Vec<T> {
    let mut vec = Vec::<u8>::new();
    for byte in &bytes[offset..] { vec.push(*byte); }
    unsafe {
        let (_, values, _) = vec.align_to::<T>();
        values.to_vec()
    }
}

fn show_ints(stuff: &[i32], offset: usize) {
    for _ in 0..offset { print!("   "); }
    for &x in stuff {
        if -100000 < x && x < 100000 && x != 0
        { print!("{x:+011} "); }
        else { print!("----------- ") };
    }
    println!("");
}

fn show_reals(stuff: &[f64], offset: usize) {
    for _ in 0..offset { print!("   "); }
    for &x in stuff {
        if x.abs() < 1.0e5 && x.abs() > 1.0e-5 && x != 0.0
        { print!("{x:+023.2e} "); }
        else { print!("----------------------- "); };
    }
    println!("");
}

fn show_reals_32(stuff: &[f32], offset: usize) {
    for _ in 0..offset { print!("   "); }
    for &x in stuff {
        if x.abs() < 1.0e5 && x.abs() > 1.0e-5 && x != 0.0
        { print!("{x:+011.2e} "); }
        else { print!("----------- "); };
    }
    println!("");
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    match args.history_file_type {
        HistoryFileType::Standard => standard(args),
        HistoryFileType::Custom   => custom(args),
    }
}
