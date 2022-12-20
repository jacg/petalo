use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::SeekFrom;

use binrw::{binrw, io::Seek};
use binrw::BinReaderExt;

use super::Point192;

#[binrw]
#[derive(Debug)]
pub(crate) struct Record {
    pub n_pairs: u8,

    #[br(count = n_pairs)]
    pub pairs: Vec<Photon>,
}

pub struct Event {
    pub blues: Vec<Photon>,
    pub pinks: Vec<Photon>,
}

impl Event {
    pub fn read(file: &mut File) -> Result<Self, Box<dyn Error>> {
        let Record { pairs: blues, .. } = file.read_le::<Record>()?;
        let Record { pairs: pinks, .. } = file.read_le::<Record>()?;
        Ok(Event { blues, pinks })
    }
}

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

pub fn iterate_events(file: &mut File) -> Result<impl Iterator<Item = Event> + '_, Box<dyn Error>> {
    crate::skip_header(file)?;
    Ok(std::iter::from_fn(|| Event::read(file).ok()))
}

// Does not use `iterate_events` as here we might want to inspect details that it ignores
pub fn show_file(file: impl AsRef<Path>, stop_after: Option<usize>) -> Result<(), Box<dyn Error>> {
    let mut file = File::open(file)?;
    file.seek(SeekFrom::Start(2_u64.pow(15)))?;
    let mut count = 0;
    let mut blue = true;
    while let Ok(Record { n_pairs: n_photons, pairs }) = file.read_le::<Record>() {
        if blue { println!("============================================================"); }
        if let Some(stop) = stop_after { if count >= stop { break } }; count += 1;
        println!("------ N {} photons: {n_photons} --------", if blue { "blue" } else { "pink" });
        for Photon { pos: Point192 { x, y, z },  energy, travel_distance, scatters_in_object: so, scatters_in_collimator: sc, weight } in pairs {
            let t = travel_distance / 0.03;
            println!("({x:7.2} {y:7.2} {z:7.2})   E:{energy:7.2}   t:{t:4.1} ps  w: {weight:4.2}   scatters obj:{so:2} col:{sc:2}");
        }
        blue = ! blue;
    }
    Ok(())
}
