use binrw::binrw;
use binrw::io::{Cursor, Seek, SeekFrom};
use binrw::{BinReaderExt, BinWriterExt};

#[binrw]
#[derive(Debug)]
struct Vec3 {
    x: f32,
    y: f32,
    z: f32,
}

#[binrw]
#[derive(Debug)]
struct Mesh {
    #[bw(calc = (vertices.len() - 1) as u32)]
    vertex_count: u32,

    #[br(count = vertex_count + 1)]
    vertices: Vec<Vec3>,
}

fn main() {
    let mesh = Mesh {
        vertices: vec![ Vec3 { x:  1.0, y:  2.0, z:  3.0, },
                        Vec3 { x: 11.0, y: 12.0, z: 13.0, },
                        Vec3 { x: 21.0, y: 22.0, z: 23.0, }, ],
    };

    let mut writer = Cursor::new(Vec::new());
    writer.write_le(&mesh).unwrap();

    let mut reader = writer;
    reader.seek(SeekFrom::Start(0)).unwrap();

    let _: Mesh = dbg!(reader.read_le().unwrap());
}
