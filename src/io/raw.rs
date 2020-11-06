/// Read / write float arrays as raw binary

use std::fs::File;
use std::io::{Write, Read};

pub fn write<'a>(data: impl Iterator<Item = &'a f32>, path: &std::path::PathBuf) -> std::io::Result<()> {
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
        write(original_data.iter(), &file_path)?;

        // Read data back from file
        let reloaded_data: Vec<_> = read_bin(&file_path)?
            .collect::<Result<_, _>>()?;

        // Check that roundtrip didn't corrupt the data
        assert_eq!(original_data, reloaded_data);
        Ok(())
    }
}