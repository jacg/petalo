use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::io::SeekFrom;

use binrw::{binrw, io::Seek};

pub mod standard;
pub mod custom;


#[derive(Debug, PartialEq, Eq)]
pub enum PhotonColour { Blue, Pink }


#[binrw]
#[derive(Debug)]
pub struct Point192 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
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
