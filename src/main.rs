use std::error::Error;
use std::io;
use std::process;

use csv;


use serde::Deserialize;
type F = f64;

use petalo::weights as pet;

#[derive(Debug, Deserialize)]
struct Event {
    event_id: usize,
    true_r1: F, true_phi1: F, true_z1: F, true_t1: F,
    true_r2: F, true_phi2: F, true_z2: F, true_t2: F,
    reco_r1: F, reco_phi1: F, reco_z1: F,
    reco_r2: F, reco_phi2: F, reco_z2: F,
}

type VEC = Vec<F>;

struct Events {
    x1: VEC, y1: VEC, z1: VEC, t1: VEC,
    x2: VEC, y2: VEC, z2: VEC, t2: VEC,
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

    pet::LOR::new(t1, t2,
                  pet::Point::new(x1, y1, z1),
                  pet::Point::new(x2, y2, z2),
    )

}


fn main() -> Result<(), Box<dyn Error>> {
    let mut measured_lors = vec![];

    // Build the CSV reader and iterate over each record.
    let file = std::fs::File::open("run_fastfastmc_1M_events.csv")?;
    let buf_reader = std::io::BufReader::new(file);
    let mut rdr = csv::Reader::from_reader(buf_reader);
    println!("Starting to read");
    for result in rdr.deserialize() {
        // The iterator yields Result<StringRecord, Error>, so we check the
        // error here..
        let event: Event = result?;
        measured_lors.push(event_to_lor(event));
    }
    println!("Collected events");
    let vbox = pet::VoxelBox::new((180.0, 180.0, 180.0), (60, 60, 60));
    pet::Image::mlem(vbox, &measured_lors);
    Ok(())
}


// use petalo::weights::{Point, VoxelBox};
// use petalo::visualize;

// fn main() {

//     let p1 = Point::new(-265.1371093069, 76.0, -200.0);
//     let p2 = Point::new( 276.002,         0.0,  200.0);
//     let vbox = VoxelBox::new((120.0, 100.0, 80.0), (50, 40, 30));

//     visualize::lor_weights(10.00, 10.03, p1, p2, vbox);
// }
