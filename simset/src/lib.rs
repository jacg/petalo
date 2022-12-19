use std::error::Error;
use std::fs::File;
use std::path::Path;

use binrw::binrw;
use binrw::io::{Seek, SeekFrom};
use binrw::BinReaderExt;

pub mod standard;

#[binrw]
#[derive(Debug)]
struct CustomRecord {
    n_pairs: u8,

    #[br(count = n_pairs)]
    pairs: Vec<custom::Photon>,
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
    while let Ok(CustomRecord { n_pairs: n_photons, pairs }) = file.read_le::<CustomRecord>() {
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

pub fn compare(file1: &Path, file2: &Path, stop_after: Option<usize>) -> Result<(), Box<dyn Error>> {
    let mut file1 = File::open(file1)?; file1.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut file2 = File::open(file2)?; file2.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    let mut read1 = || standard::Event::read(&mut file1);
    let mut read2 = || standard::Event::read(&mut file2);
    type D = standard::Event;
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
