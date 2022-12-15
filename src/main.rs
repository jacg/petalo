use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::path::{Path, PathBuf};

// ----- CLI --------------------------------------------------
use clap::Parser;

/// Try to parse SimSET (custom) history files
#[derive(Parser, Debug)]
struct Args {
    /// History file to parse
    file: PathBuf,

    /// Second history file for record-by-record comparison
    #[arg(short = 'c', long)]
    compare: Option<PathBuf>,

    /// Maximum number of records to parse
    #[arg(short = 'n', long)]
    stop_after: Option<usize>,

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
struct Decay {
    pos       : Point192, // 24
    weight    : f64,      //  8
    time      : f64,      //  8
    decay_type: u64,      //  8  // Needs to be mapped to an enum (PhgEn_DecayTypeTy)
} // bytes wasted: probably all 48 of them

#[binrw]
#[derive(Debug)]
struct Photon { // Both blue and pink?                        -- comments from struct definicion in SimSET: Photon.h --
    location              : Point96,  // 123456789aba  12  12 photon current location or, in detector list mode, detection position
    angle                 : Point96,  // 123456789aba  12  24 photon durrent direction.  perhaps undefined in detector list mode.
    flags                 : u8,       // 1              1  25 Photon flags
    #[br(pad_before = 7)]             //  2345678       7  32
    weight                : f64,      // 12345678       8  40 Photon's weight
    energy                : f32,      // 1234           4  44 Photon's energy
    #[br(pad_before = 4)]             //     5678       4  48
    time                  : f64,      // 12345678       3  56 In seconds
    transaxial_position   : f32,      // 1234           4  60 For SPECT, transaxial position on "back" of collimator/detector
    azimuthal_angle_index : u16,      // 12             2  62 For SPECT/DHCI, index of collimator/detector angle
    #[br(pad_before = 2)]             //   34           2  64
    detector_angle        : f32,      // 1234           4  68 For SPECT, DHCI, angle of detector
    detector_crystal      : u32,      // 1234           4  72 for block detectors, the crystal number for detection
} // bytes wasted: 17 on padding, 10 on SPECT, 12 on undefined-in-LM direction, 4 on block detectors; 43/72 = 60%

#[binrw]
#[derive(Debug)]
enum Standard {

    #[br(magic = 1_u8)]
    Decay(Decay),

    // according to struct PHG_DetectedPhoton in file Photon.h
    #[br(magic = 2_u8)]
    Photon(Photon),
}

#[binrw]
#[derive(Debug)]
struct Custom {
    n_pairs: u8,

    #[br(count = n_pairs)]
    pairs: Vec<Entry>,
}

#[binrw]
#[derive(Debug)]
struct Entry {
    pos: Point192,
    // x-cos
    // y-cos
    // z-cos
    // scatters in object
    // scatters in collimator
    // decay_weight: f64,
    // weight: f64,
    energy: f64,
    travel_distance: f64,
    // decay_pos: Point192,
    // decay_time: f64,
    // decay_type: u32,

    // LOR
    // transaxial_distance: f64
    // azimuthal angle index
    // axial_position: f64,

    // detector_pos: Point192,
    // detector_angle: f64,
    // detector crystal

    // num_detector_interactions: i32, // Supposed to be 1 char, but it looks like alignment sometimes pads it!
    // num_detector_interactions_written: i32, // ? Will be zero if detector_interactions not stored

    // #[br(count = num_detector_interactions_written)]
    // detector_interactions: Vec<DetectorInteraction>,

}

#[binrw]
#[derive(Debug)]
struct Point192 {
    x: f64,
    y: f64,
    z: f64,
}

#[binrw]
#[derive(Debug)]
struct Point96 {
    x: f32,
    y: f32,
    z: f32,
}

#[binrw]
#[derive(Debug)]
struct DetectorInteraction {
    pos: Point192,
    energy_deposited: f64,
    is_active: u8,
    padding: u8, // Is this here because of alignment before writing?
}

fn custom(args: Args) -> Result<(), Box<dyn Error>> {
    let mut file = File::open(args.file)?;
    file.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    while let Ok(Custom { n_pairs, pairs }) = file.read_le::<Custom>() {
        if let Some(stop) = args.stop_after { if count >= stop { break } }; count += 1;
        println!("------ N pairs: {n_pairs} --------");
        for Entry { pos: Point192 { x, y, z },  energy, travel_distance } in pairs {
            println!("({x:7.2} {y:7.2} {z:7.2})   E:{energy:7.2}   dist:{travel_distance:6.2} cm");
        }
    }
    Ok(())
}

fn standard(args: Args) -> Result<(), Box<dyn Error>> {
    let mut file = File::open(args.file)?;
    file.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    let mut ts = [None, None];
    while let Ok(record) = file.read_le::<Standard>() {
        match record {
            Standard::Decay(Decay { pos: Point192 { x, y, z }, weight, time, decay_type }) => {
                if let Some(stop) = args.stop_after { if count >= stop { break } }; count += 1;
                ts = [None, None];
                println!("\n===================================================================================");
                println!("({x:6.2} {y:6.2} {z:6.2})    t: {time:5.2}  s   w:{weight:4.1}   type:{decay_type}");
            },
            Standard::Photon(Photon { location: Point96 { x, y, z }, flags, weight, energy, time, .. }) => {
                let time = time * 1e12;
                ts[blue_or_pink_index(flags)] = Some(time);
                let (c, s, t) = interpret_flags(flags);
                println!("({x:6.2} {y:6.2} {z:6.2})    + {time:6.1} ps   w:{weight:4.1}  E: {energy:6.2}   {c} {s} {t} {flags:02}");
                if let [Some(t1), Some(t2)] = ts {
                    let dtof = t1 - t2;
                    let dx = dtof * 0.3 /* mm / ps */ / 2.0;
                    println!("dTOF: {dtof:6.1} ps    ->    dx {dx:6.2} mm");
                }
            },
        }
    }
    println!("Stopping after {count} records");
    Ok(())
}

fn compare(file1: &Path, file2: &Path, args: &Args) -> Result<(), Box<dyn Error>> {
    let mut file1 = File::open(file1)?; file1.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut file2 = File::open(file2)?; file2.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    let mut read1 = || StandardDecayWithPhotons::read(&mut file1);
    let mut read2 = || StandardDecayWithPhotons::read(&mut file2);
    type D = StandardDecayWithPhotons;
    while let (Ok(D { decay: ldecay, photons: lphotons }), Ok(D { decay: rdecay, photons: rphotons })) = (read1(), read2()) {
        println!();
        if let Some(stop) = args.stop_after {
            if count >= stop { break };
            count += 1;
        }
        for (Photon { energy: el, .. }, Photon { energy: er, .. }) in lphotons.into_iter().zip(rphotons) {
            println!("{el:7.2} - {er:7.2} = {:7.2}", el - er);
        }
    }
    Ok(())
}

fn blue_or_pink_index(flag: u8) -> usize { (flag & 0x1) as usize }
fn interpret_flags(flag: u8) -> (&'static str, &'static str, &'static str) {
    let c = if (flag & 0x1) == 0x1 { "blue"    } else { "pink"    };
    let s = if (flag & 0x2) == 0x2 { "scatter" } else { "-------" };
    let t = if (flag & 0x4) == 0x4 { "primary" } else { "-------" };
    (c, s, t)
}

struct StandardDecayWithPhotons {
    decay: Decay,
    photons: Vec<Photon>,
}

const DECAY_SIZE_INCLUDING_MAGIC_: i64 = (std::mem::size_of::<Decay>() + 1) as i64;

impl StandardDecayWithPhotons {
    fn read(file: &mut File) -> Result<Self, Box<dyn Error>> {
        use Standard::*;
        let mut photons = vec![];
        let       Decay (decay ) = file.read_le::<Standard>()? else { panic!("Expected decay") };
        while let Photon(photon) = file.read_le::<Standard>()? { photons.push(photon); }
        file.seek(std::io::SeekFrom::Current(-DECAY_SIZE_INCLUDING_MAGIC_))?;
        Ok(Self { decay, photons })
    }
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

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if let Some(ref file2) = args.compare { compare(&args.file, &file2, &args) }
    else {
        match args.history_file_type {
            HistoryFileType::Standard => standard(args),
            HistoryFileType::Custom   => custom(args),
        }
    }
}
