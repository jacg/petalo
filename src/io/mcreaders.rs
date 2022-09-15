use std::path::Path;
use crate::io::hdf5::read_table;
use crate::Time;

use crate::config::mlem::Bounds;

use geometry::units::ns;

// Use otherwise pointless module to allow nonstandard_style in constants
// generated by hdf5 derive macro
pub use grr::*;
#[allow(nonstandard_style)]
mod grr {

    #[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
    #[repr(C)]
    pub struct MCVertex {
        pub event_id: u32,
        pub track_id: u32,
        pub parent_id: u32,
        pub x: f32,
        pub y: f32,
        pub z: f32,
        pub t: f32,
        pub moved: f32,
        pub pre_KE: f32,
        pub post_KE: f32,
        pub deposited: u32,
        pub process_id: u32, // NB these may differ across
        pub volume_id: u32,  // different files
    }

    #[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
    #[repr(C)]
    pub struct MCSensorHit {
        pub event_id: u32,
        pub sensor_id: u32,
        pub time: f32,
    }

    #[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
    #[repr(C)]
    pub struct MCQtot {
        pub event_id: u32,
        pub sensor_id: u32,
        pub charge: u32,
    }
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct MCPrimary {
    pub event_id: u32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub vx: f32,
    pub vy: f32,
    pub vz: f32,
}

#[derive(hdf5::H5Type, Clone, PartialEq, Debug)]
#[repr(C)]
pub struct MCSensorXYZ {
    pub sensor_id: u32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

// Combines individual SensorHit times with total charge.
#[derive(Clone, Debug)]
pub struct MCQT {
    pub event_id: u32,
    pub sensor_id: u32,
    pub q: u32,
    pub t: Time,
}

pub fn read_vertices(filename: &Path, range: Bounds<usize>) -> hdf5::Result<Vec<MCVertex>> {
    Ok(array_to_vec(read_table::<MCVertex>(&filename, "MC/vertices", range)?))
}

pub fn read_sensor_hits(filename: &Path, range: Bounds<usize>) -> hdf5::Result<Vec<MCSensorHit>> {
    Ok(array_to_vec(read_table::<MCSensorHit>(&filename, "MC/waveform", range)?))
}

pub fn read_primaries(filename: &Path, range: Bounds<usize>) -> hdf5::Result<Vec<MCPrimary>> {
    Ok(array_to_vec(read_table::<MCPrimary>(&filename, "MC/primaries", range)?))
}

pub fn read_sensor_xyz(filename: &Path) -> hdf5::Result<Vec<MCSensorXYZ>> {
    Ok(array_to_vec(read_table::<MCSensorXYZ>(&filename, "MC/sensor_xyz", Bounds::none())?))
}

pub fn read_qts(infile: &Path, range: Bounds<usize>) -> hdf5::Result<Vec<MCQT>> {
    // Read charges and waveforms
    let qs = read_table::<MCQtot     >(&infile, "MC/total_charge", range.clone())?;
    let ts = read_table::<MCSensorHit>(&infile, "MC/waveform"    , range        )?;
    Ok(combine_tables(qs, ts))
}

// TODO Is there really no simpler way?
pub fn array_to_vec<T: Clone>(array: ndarray::Array1<T>) -> Vec<T> {
    let mut vec = vec![];
    vec.extend_from_slice(array.as_slice().unwrap());
    vec
}

fn combine_tables(qs: ndarray::Array1<MCQtot>, ts: ndarray::Array1<MCSensorHit>) -> Vec<MCQT> {
    let mut qts = vec![];
    let mut titer = ts.iter();
    for &MCQtot{ event_id, sensor_id, charge:q} in qs.iter() {
        for &MCSensorHit{ event_id: te, sensor_id: ts, time:t} in titer.by_ref() {
            if event_id == te && sensor_id == ts {
                qts.push(MCQT{ event_id, sensor_id, q, t: ns(t) });
                break;
            }
        }
    }
    qts
}

// -------- Tests ------

#[cfg(test)]
mod test_mcreaders {
    use super::*;

    fn get_test_file() -> &'static Path {
        Path::new("src/io/mcead_test.h5")
    }

    #[test]
    fn test_read_vertices() {
        let test_file = get_test_file();

        let vertices = read_vertices(test_file, Bounds::none());
        let track_ids: Vec<u32> = vertices.unwrap().into_iter().map(|vrt| vrt.track_id).collect();
        assert_eq!(track_ids, [1, 1, 2, 4, 1]);
    }

    #[test]
    fn test_read_sensor_hits() {
        let test_file = get_test_file();

        let sensor_hits = read_sensor_hits(test_file, Bounds::none());
        assert_eq!(sensor_hits.unwrap().len(), 8098);
    }

    #[test]
    fn test_read_primaries() {
        let test_file = get_test_file();

        let primaries = read_primaries(test_file, Bounds::none());
        assert_eq!(primaries.unwrap().len(), 5);
    }

    #[test]
    fn test_read_sensor_xyz() {
        let test_file = get_test_file();

        let sensor_xyz = read_sensor_xyz(test_file);
        assert_eq!(sensor_xyz.unwrap().len(), 47428);
    }
}