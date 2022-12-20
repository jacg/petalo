use std::error::Error;
use std::fs::File;
use std::path::Path;

use binrw::{binrw, io::Seek};
use binrw::{BinReaderExt, BinRead};

use units::{cm, Length, Time};

use crate::{f64_as_s as s, skip_header};

#[derive(Debug, BinRead)]
pub(crate) enum Record {
    #[br(magic = 1_u8)] Decay (Decay),
    #[br(magic = 2_u8)] Photon(Photon),
    // according to struct PHG_DetectedPhoton in file Photon.h
}

pub struct Event {
    pub decay: Decay,
    pub photons: Vec<Photon>,
}

const DECAY_SIZE_INCLUDING_MAGIC_: i64 = (std::mem::size_of::<Decay>() + 1) as i64;

impl Event {
    pub fn read(file: &mut File) -> Result<Self, Box<dyn Error>> {
        let mut photons = vec![];
        let       Record::Decay (decay ) = file.read_le::<Record>()? else { panic!("Expected decay") };
        while let Record::Photon(photon) = file.read_le::<Record>()? { photons.push(photon); }
        file.seek(std::io::SeekFrom::Current(-DECAY_SIZE_INCLUDING_MAGIC_))?;
        Ok(Self { decay, photons })
    }
}

#[binrw]
#[derive(Debug)]
pub struct Decay {
    pub x         : f64,      //  8
    pub y         : f64,      //  8
    pub z         : f64,      //  8
    pub weight    : f64,      //  8
    pub time      : f64,      //  8
    pub decay_type: u64,      //  8  // Needs to be mapped to an enum (PhgEn_DecayTypeTy)
} // bytes wasted: probably all 48 of them

#[derive(Debug, BinRead)]
pub struct Photon { // Both blue and pink?                        -- comments from struct definicion in SimSET: Photon.h --
    #[br(map = cm)] pub x     : Length,   // 1234      4   4 photon current location or, in detector list mode, detection position
    #[br(map = cm)] pub y     : Length,   // 1234      4   8 photon current location or, in detector list mode, detection position
    #[br(map = cm)] pub z     : Length,   // 1234      4  12 photon current location or, in detector list mode, detection position
    pub angle_x               : f32,      // 1234      4  16 photon durrent direction.  perhaps undefined in detector list mode.
    pub angle_y               : f32,      // 1234      4  20 photon durrent direction.  perhaps undefined in detector list mode.
    pub angle_z               : f32,      // 1234      4  24 photon durrent direction.  perhaps undefined in detector list mode.
    pub flags                 : u8,       // 1         1  25 Photon flags
    #[br(pad_before = 7)]                 //  2345678  7  32
    pub weight                : f64,      // 12345678  8  40 Photon's weight
    pub energy                : f32,      // 1234      4  44 Photon's energy
    #[br(pad_before = 4)]                 //     5678  4  48
    #[br(map = s)] pub time   : Time,     // 12345678  3  56 In seconds
    pub transaxial_position   : f32,      // 1234      4  60 For SPECT, transaxial position on "back" of collimator/detector
    pub azimuthal_angle_index : u16,      // 12        2  62 For SPECT/DHCI, index of collimator/detector angle
    #[br(pad_before = 2)]                 //   34      2  64
    pub detector_angle        : f32,      // 1234      4  68 For SPECT, DHCI, angle of detector
    pub detector_crystal      : u32,      // 1234      4  72 for block detectors, the crystal number for detection
} // bytes wasted: 17 on padding, 10 on SPECT, 12 on undefined-in-LM direction, 4 on block detectors; 43/72 = 60%

impl Photon {
    pub fn colour(&self) -> super::PhotonColour {
        if self.flags & 0x1 == 0x1 { super::PhotonColour::Blue }
        else                       { super::PhotonColour::Pink }
    }
}

pub fn iterate_events(file: &mut File) -> Result<impl Iterator<Item = Event> + '_, Box<dyn Error>> {
    skip_header(file)?;
    Ok(std::iter::from_fn(|| Event::read(file).ok()))
}

// Does not use `iterate_events` as here we might want to inspect details that it ignores
pub fn show_file(file: impl AsRef<Path>, stop_after: Option<usize>) -> Result<(), Box<dyn Error>> {
    let mut file = File::open(file)?;
    super::skip_header(&mut file)?;
    let mut count = 0;
    let mut _ts = [None, None];
    while let Ok(decay) = Event::read(&mut file) {// file.read_le::<Standard>() {

        // ----- Process the decay data -------------------------------------------
        let Decay { x, y, z, weight, time, decay_type } = decay.decay;
        if let Some(stop) = stop_after { if count >= stop { break } }; count += 1;
        _ts = [None, None];
        println!("=====================================================================================");
        println!("({x:8.2} {y:8.2} {z:8.2})    t: {time:5.2}  s   w:{weight:6.2}   type:{decay_type}");

        // ----- Process each photon associated with the decay --------------------
        let mut seen_pink = false;
        for Photon { x, y, z, flags, weight, energy, time, .. } in decay.photons {
            _ts[blue_or_pink_index(flags)] = Some(time);
            let (c, s, t) = interpret_flags(flags);
            if c == "pink" && !seen_pink {
                seen_pink = true;
                println!();
            }
            println!("({x:7.2?} {y:7.2?} {z:7.2?})    + {time:6.1?}   w:{weight:6.2}  E: {energy:6.2}   {c} {s} {t} {flags:02}");
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

fn blue_or_pink_index(flag: u8) -> usize { (flag & 0x1) as usize }

fn interpret_flags(flag: u8) -> (&'static str, &'static str, &'static str) {
    let c = if (flag & 0x1) == 0x1 { "blue"    } else { "pink"    };
    let s = if (flag & 0x2) == 0x2 { "scatter" } else { "-------" };
    let t = if (flag & 0x4) == 0x4 { "primary" } else { "-------" };
    (c, s, t)
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
