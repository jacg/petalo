use super::{Point192, binrw};

#[binrw]
#[derive(Debug)]
pub(crate) struct Record {
    pub n_pairs: u8,

    #[br(count = n_pairs)]
    pub pairs: Vec<Photon>,
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
