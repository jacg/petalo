/// Read / write float arrays as raw binary

use std::fs::File;
use std::io::{Write, Read, BufWriter};

pub fn write(data: impl Iterator<Item = f32>, path: &std::path::Path) -> std::io::Result<()> {
    let file = File::create(path)?;
    let mut buf = BufWriter::new(file);
    for datum in data {
        let bytes = datum.to_le_bytes();
        // TODO: Clippy suggests write -> write_all, but that slows down the
        // writer by a factor of 100. However, like this, we're not checking
        // whether the whole buffer was written.
        buf.write(&bytes)?;
    }
    Ok(())
}

type IORes<T> = std::io::Result<T>;
// TODO: reading is a factor of 60 slower than writing
pub fn read<'a>(path: &std::path::Path) -> IORes<impl Iterator<Item = IORes<f32>> + 'a> {
    let mut file = File::open(path)?;
    let mut buffer = [0; 4];

    Ok(std::iter::from_fn(move || {
        use std::io::ErrorKind::UnexpectedEof;
        match file.read_exact(&mut buffer) {
            Ok(()) => Some(Ok(f32::from_le_bytes(buffer))),
            Err(e) if e.kind() == UnexpectedEof => None,
            Err(e) => Some(Err(e)),
        }
    }))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn raw_io_roundtrip() -> std::io::Result<()> {
        use tempfile::tempdir;
        #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};

        // Harmless temporary location for output file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.bin");

        // Some test data
        let original_data = vec![1.23, 4.56, 7.89];

        // Write data to file
        write(original_data.iter().copied(), &file_path)?;

        // Read data back from file
        let reloaded_data: Vec<_> = read(&file_path)?
            .collect::<Result<_, _>>()?;

        // Check that roundtrip didn't corrupt the data
        assert_eq!(original_data, reloaded_data);
        Ok(())
    }
}
