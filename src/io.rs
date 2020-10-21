use std::error::Error;
use serde::{Serialize, Deserialize};

use crate::weights as pet;

type F = pet::Length;


#[derive(Debug, Serialize, Deserialize, Clone)]
struct Event {
    event_id: usize,
    true_r1: F, true_phi1: F, true_z1: F, true_t1: F,
    true_r2: F, true_phi2: F, true_z2: F, true_t2: F,
    reco_r1: F, reco_phi1: F, reco_z1: F,
    reco_r2: F, reco_phi2: F, reco_z2: F,
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
