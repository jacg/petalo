use std::error::Error;
use std::fs::File;
use std::path::Path;

use binrw::binrw;
use binrw::io::{Seek, SeekFrom};
use binrw::BinReaderExt;


#[binrw]
#[derive(Debug)]
enum Standard {
    #[br(magic = 1_u8)] Decay (standard::Decay),
    #[br(magic = 2_u8)] Photon(standard::Photon),
    // according to struct PHG_DetectedPhoton in file Photon.h
}

#[binrw]
#[derive(Debug)]
struct Custom {
    n_pairs: u8,

    #[br(count = n_pairs)]
    pairs: Vec<custom::Photon>,
}

mod standard {
    use super::{Point96, Point192, binrw};

    #[binrw]
    #[derive(Debug)]
    pub struct Decay {
        pub pos       : Point192, // 24
        pub weight    : f64,      //  8
        pub time      : f64,      //  8
        pub decay_type: u64,      //  8  // Needs to be mapped to an enum (PhgEn_DecayTypeTy)
    } // bytes wasted: probably all 48 of them

    #[binrw]
    #[derive(Debug)]
    pub struct Photon { // Both blue and pink?                        -- comments from struct definicion in SimSET: Photon.h --
        pub location              : Point96,  // 123456789aba  12  12 photon current location or, in detector list mode, detection position
        pub angle                 : Point96,  // 123456789aba  12  24 photon durrent direction.  perhaps undefined in detector list mode.
        pub flags                 : u8,       // 1              1  25 Photon flags
        #[br(pad_before = 7)]                 //  2345678       7  32
        pub weight                : f64,      // 12345678       8  40 Photon's weight
        pub energy                : f32,      // 1234           4  44 Photon's energy
        #[br(pad_before = 4)]                 //     5678       4  48
        pub time                  : f64,      // 12345678       3  56 In seconds
        pub transaxial_position   : f32,      // 1234           4  60 For SPECT, transaxial position on "back" of collimator/detector
        pub azimuthal_angle_index : u16,      // 12             2  62 For SPECT/DHCI, index of collimator/detector angle
        #[br(pad_before = 2)]                 //   34           2  64
        pub detector_angle        : f32,      // 1234           4  68 For SPECT, DHCI, angle of detector
        pub detector_crystal      : u32,      // 1234           4  72 for block detectors, the crystal number for detection
    } // bytes wasted: 17 on padding, 10 on SPECT, 12 on undefined-in-LM direction, 4 on block detectors; 43/72 = 60%

    impl Photon {
        pub fn colour(&self) -> super::PhotonColour {
            if self.flags & 0x1 == 0x1 { super::PhotonColour::Blue }
            else                       { super::PhotonColour::Pink }
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum PhotonColour { Blue, Pink }

mod custom {
    use super::{Point192, binrw};

    #[binrw]
    #[derive(Debug)]
    pub struct Photon {
        pub pos: Point192,
        // x-cos
        // y-cos
        // z-cos
        pub scatters_in_object: u32,
        pub scatters_in_collimator: u32,
        // decay_weight: f64,
        pub weight: f64,
        pub energy: f64,
        pub travel_distance: f64,
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
}
#[binrw]
#[derive(Debug)]
pub struct Point192 {
    x: f64,
    y: f64,
    z: f64,
}

#[binrw]
#[derive(Debug)]
pub struct Point96 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[binrw]
#[derive(Debug)]
struct DetectorInteraction {
    pos: Point192,
    energy_deposited: f64,
    is_active: u8,
    padding: u8, // Is this here because of alignment before writing?
}

pub fn custom(file: impl AsRef<Path>, stop_after: Option<usize>) -> Result<(), Box<dyn Error>> {
    let mut file = File::open(file)?;
    file.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    let mut blue = true;
    while let Ok(Custom { n_pairs: n_photons, pairs }) = file.read_le::<Custom>() {
        if blue { println!("============================================================"); }
        if let Some(stop) = stop_after { if count >= stop { break } }; count += 1;
        println!("------ N {} photons: {n_photons} --------", if blue { "blue" } else { "pink" });
        for custom::Photon { pos: Point192 { x, y, z },  energy, travel_distance, scatters_in_object: so, scatters_in_collimator: sc, weight } in pairs {
            let t = travel_distance / 0.03;
            println!("({x:7.2} {y:7.2} {z:7.2})   E:{energy:7.2}   t:{t:4.1} ps  w: {weight:4.2}   scatters obj:{so:2} col:{sc:2}");
        }
        blue = ! blue;
    }
    Ok(())
}

pub fn skip_header(file: &mut File) -> std::io::Result<u64> {
    file.seek(SeekFrom::Start(2_u64.pow(15)))
}

pub fn standard(file: impl AsRef<Path>, stop_after: Option<usize>) -> Result<(), Box<dyn Error>> {
    use standard as st;
    let mut file = File::open(file)?;
    skip_header(&mut file)?;
    let mut count = 0;
    let mut _ts = [None, None];
    while let Ok(decay) = StandardDecayWithPhotons::read(&mut file) {// file.read_le::<Standard>() {

        // ----- Process the decay data -------------------------------------------
        let st::Decay { pos: Point192 { x, y, z }, weight, time, decay_type } = decay.decay;
        if let Some(stop) = stop_after { if count >= stop { break } }; count += 1;
        _ts = [None, None];
        println!("=====================================================================================");
        println!("({x:6.2} {y:6.2} {z:6.2})    t: {time:5.2}  s   w:{weight:6.2}   type:{decay_type}");

        // ----- Process each photon associated with the decay --------------------
        let mut seen_pink = false;
        for st::Photon { location: Point96 { x, y, z }, flags, weight, energy, time, .. } in decay.photons {
            let time = time * 1e12;
            _ts[blue_or_pink_index(flags)] = Some(time);
            let (c, s, t) = interpret_flags(flags);
            if c == "pink" && !seen_pink {
                seen_pink = true;
                println!();
            }
            println!("({x:6.2} {y:6.2} {z:6.2})    + {time:6.1} ps   w:{weight:6.2}  E: {energy:6.2}   {c} {s} {t} {flags:02}");
            if let [Some(t1), Some(t2)] = _ts {
                let dtof = t1 - t2;
                let _dx = dtof * 0.3 /* mm / ps */ / 2.0;
                //println!("dTOF: {dtof:6.1} ps    ->    dx {dx:6.2} mm");
            }
        }
    }
    println!("Stopping after {count} records");
    Ok(())
}

pub fn compare(file1: &Path, file2: &Path, stop_after: Option<usize>) -> Result<(), Box<dyn Error>> {
    let mut file1 = File::open(file1)?; file1.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut file2 = File::open(file2)?; file2.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    let mut read1 = || StandardDecayWithPhotons::read(&mut file1);
    let mut read2 = || StandardDecayWithPhotons::read(&mut file2);
    type D = StandardDecayWithPhotons;
    while let (Ok(D { decay: _ldecay, photons: lphotons }), Ok(D { decay: _rdecay, photons: rphotons })) = (read1(), read2()) {
        //println!();
        if let Some(stop) = stop_after {
            if count >= stop { break };
            count += 1;
        }
        for (standard::Photon { energy: el, .. }, standard::Photon { energy: er, .. }) in lphotons.into_iter().zip(rphotons) {
            //println!("{el:7.2} - {er:7.2} = {:7.2}", el - er);
            println!("{:7.2}", el - er);
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
    decay: standard::Decay,
    photons: Vec<standard::Photon>,
}

const DECAY_SIZE_INCLUDING_MAGIC_: i64 = (std::mem::size_of::<standard::Decay>() + 1) as i64;

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
