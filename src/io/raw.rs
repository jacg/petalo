/// Read / write float arrays as raw binary

use std::fs::File;
use std::io::{Write, Read, BufWriter, BufReader};

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
pub fn read<'a>(path: &std::path::Path) -> IORes<impl Iterator<Item = IORes<f32>> + 'a> {
    let file = File::open(path)?;
    let mut buf = BufReader::new(file);
    let mut buffer = [0; 4];

    Ok(std::iter::from_fn(move || {
        use std::io::ErrorKind::UnexpectedEof;
        match buf.read_exact(&mut buffer) {
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

// ----- Raw 3d image with matrix/physical size metadata --------------------------------

use binrw::{BinRead, BinReaderExt};
use binwrite::BinWrite;

#[derive(BinRead, BinWrite, PartialEq, Debug)]
// #[br(magic = b"IMG3D")] // TODO: magic not supported by writers?
#[br(big)]
#[binwrite(big)]
pub struct Image3D {
    pub pixels: [u16; 3],
    pub mm: [f32; 3],
    #[br(count = pixels[0] as usize * pixels[1] as usize * pixels[2] as usize)]
    pub data: Vec<f32>,
}

impl Image3D {
    pub fn write_to_file(&self, path: impl AsRef<std::path::Path>) -> std::io::Result<()> {
        let file = File::create(path)?;
        let mut buf = BufWriter::new(file);
        self.write(&mut buf)
    }

    pub fn read_from_file(path: impl AsRef<std::path::Path>) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut buffered = BufReader::new(file);
        Ok(buffered.read_ne().unwrap())
    }
}

// TODO: Ideally mlem::Image and Image3D would be a single type.
use super::super::mlem::Image as MLEMImage;
impl From<&MLEMImage> for Image3D {
    fn from(image: &MLEMImage) -> Self {
        let n = image.vbox.n;
        let pixels = [n[0] as u16, n[1] as u16, n[2] as u16];
        let l = image.vbox.half_width;
        let mm = [l[0]*2.0, l[1]*2.0, l[2]*2.0];
        let data = image.data.clone();
        Self { pixels, mm, data }
    }
}

impl From<&Image3D> for MLEMImage {
    fn from(image: &Image3D) -> Self {
        let [px, py, pz] = image.pixels;
        let n = (px as usize, py as usize, pz as usize);
        let [wx, wy, wz] = image.mm;
        let half_width = (wx, wy, wz);
        let vbox = crate::weights::VoxelBox::new(half_width, n);
        let data = image.data.clone();
        Self { vbox, data }
    }
}

#[cfg(test)]
use binrw::io::Cursor;

#[cfg(test)]
mod test_with_metadata {
    use super::*;

    // An example instance used for testing
    fn guinea_pig() -> Image3D {
        Image3D {
            pixels: [1, 2, 3],
            mm:     [2.0, 4.0, 9.0],
            data: (0..6).map(|n| n as f32).collect(),
        }
    }

    #[test]
    fn roundtrip_via_buffer() {
        // Create a test image
        let original = guinea_pig();

        // Serialize the image
        let mut bytes = vec![];
        original.write(&mut bytes).unwrap();
        println!("{:?}", original);

        // Deserialize
        let mut reader = Cursor::new(bytes);
        let recovered: Image3D = reader.read_ne().unwrap();

        // Verify
        assert_eq!(recovered, original);
    }

    #[test]
    fn roundtrip_via_file() -> std::io::Result<()> {
        use tempfile::tempdir;
        #[allow(unused)] use pretty_assertions::{assert_eq, assert_ne};

        // Harmless temporary location for output file
        let dir = tempdir()?;
        let file_path = dir.path().join("test.bin");

        // Some test data
        let original = guinea_pig();

        // Write data to file
        original.write_to_file(&file_path)?;

        // Read data back from file
        let recovered = Image3D::read_from_file(&file_path)?;

        // Check that roundtrip didn't corrupt the data
        assert_eq!(original, recovered);
        Ok(())
    }

    #[test]
    fn roundtrip_via_mlem_image() {
        let original = guinea_pig();
        let converted = crate::mlem::Image::from(&original);
        let recovered = Image3D::from(&converted);
        assert_eq!(original, recovered);
    }

}

// ----- Proofs of concept ---------------------------------------------------------------
#[cfg(test)]
mod test_br_enum {
    use super::*;

    #[test]
    fn test_br_enum() {
        let point = Point::read(&mut Cursor::new(b"\x80\x02\xe0\x01")).unwrap();
        assert_eq!(point, Point(640, 480));

        #[derive(BinRead, Debug, PartialEq)]
        #[br(little)]
        struct Point(i16, i16);

        #[derive(BinRead, Debug, PartialEq)]
        #[br(big, magic = b"SHAP")]
        enum Shape {
            #[br(magic(0u8))] Rect { left: i16, top: i16, right: i16, bottom: i16 },
            #[br(magic(1u8))] Oval { origin: Point, rx: u8, ry: u8 }
        }

        let oval = Shape::read(&mut Cursor::new(b"SHAP\x01\x80\x02\xe0\x01\x2a\x15")).unwrap();
        assert_eq!(oval, Shape::Oval { origin: Point(640, 480), rx: 42, ry: 21 });

        let rect = Shape::read(&mut Cursor::new(b"SHAP\x00\x00\x01\x00\x02\x00\x03\x00\x04")).unwrap();
        assert_eq!(rect, Shape::Rect { left: 1, top: 2, right: 3, bottom: 4 });
    }
}
