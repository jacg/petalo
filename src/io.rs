use std::error::Error;
use serde::{Serialize, Deserialize};

use crate::weights as pet;

#[derive(Debug, Serialize, Deserialize, Clone)]
struct Event {
    event_id: usize,
    true_r1: f32, true_phi1: f32, true_z1: f32, true_t1: f32,
    true_r2: f32, true_phi2: f32, true_z2: f32, true_t2: f32,
    reco_r1: f32, reco_phi1: f32, reco_z1: f32,
    reco_r2: f32, reco_phi2: f32, reco_z2: f32,
}

    fn event_to_lor(
        Event{ reco_r1: r1, reco_phi1: phi1, reco_z1: z1, true_t1: t1,
               reco_r2: r2, reco_phi2: phi2, reco_z2: z2, true_t2: t2, .. }: Event)
    -> pet::LOR
{
    let x1 = r1 * phi1.cos();
    let y1 = r1 * phi1.sin();
    let t1 = t1 * 1000.0; // turn into picoseconds

    let x2 = r2 * phi2.cos();
    let y2 = r2 * phi2.sin();
    let t2 = t2 * 1000.0; // turn into picoseconds

    pet::LOR::new(t1.into(), t2.into(),
                  pet::Point::new(x1.into(), y1.into(), z1.into()),
                  pet::Point::new(x2.into(), y2.into(), z2.into()),
    )

}

pub fn read_lors(filename: &str) -> Result<Vec<pet::LOR>, Box<dyn Error>> {
    let mut measured_lors = vec![];
    {
        let mut file = std::io::BufReader::new(std::fs::File::open(filename)?);
        let events: Vec<Event> = bincode::deserialize_from(&mut file)?;
        for event in events {
            measured_lors.push(event_to_lor(event.clone()));
        }
    }
    Ok(measured_lors)
}

// ----------- read/write float arrays as raw binary --------------------------
use std::fs::File;
use std::io::{Write, Read};

pub fn write_bin<'a>(data: impl Iterator<Item = &'a f32>, path: &std::path::PathBuf) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for datum in data {
        let bytes = datum.to_le_bytes();
        file.write(&bytes)?;
    }
    Ok(())
}

type IORes<T> = std::io::Result<T>;

pub fn read_bin<'a>(path: &std::path::PathBuf) -> IORes<impl Iterator<Item = IORes<f32>> + 'a> {
    let mut file = File::open(path)?;
    let mut buffer = [0; 4];

    Ok(std::iter::from_fn(move || {
        use std::io::ErrorKind::UnexpectedEof;
        match file.read_exact(&mut buffer) {
            Ok(()) => Some(Ok(f32::from_le_bytes(buffer))),
            Err(e) if e.kind() == UnexpectedEof => None,
            Err(e) => return Some(Err(e)),
        }
    }))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn bin_io_roundtrip() -> std::io::Result<()> {
        use tempfile::tempdir;
        #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};

        // Harmless temporary location for output file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.bin");

        // Some test data
        let original_data = vec![1.23, 4.56, 7.89];

        // Write data to file
        write_bin(original_data.iter(), &file_path)?;

        // Read data back from file
        let reloaded_data: Vec<_> = read_bin(&file_path)?
            .collect::<Result<_, _>>()?;

        // Check that roundtrip didn't corrupt the data
        assert_eq!(original_data, reloaded_data);
        Ok(())
    }
}
