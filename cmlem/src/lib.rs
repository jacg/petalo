#![allow(nonstandard_style)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

use std::ffi::CString;
use ::std::os::raw::{c_char, c_int};

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn does_not_crash() {
        let niterations: c_int = 5;
        let TOF = false;
        let TOF_resolution: f32 = 200.00;
        let FOV_XY: f32 = 180.0;
        let FOV_Z : f32 = 180.0;
        let NXY: c_int = 60;
        let NZ : c_int = 60;    let image_size = NXY as usize * NXY as usize * NZ as usize;
        let ncoinc: c_int = 2;
        let LOR_X1: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_Y1: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_Z1: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_T1: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_X2: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_Y2: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_Z2: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let LOR_T2: *mut f32 = [-100.0, -100.0].as_mut_ptr();
        let mut SENS: Vec<f32> = vec![1.0; image_size];
        let SENS: *mut f32 = SENS.as_mut_ptr();
        let outfile_prefix = CString::new("/tmp/mlemmy").unwrap();
        let outfile_prefix: *const c_char = outfile_prefix.as_ptr();
        let out_niter: c_int = 1;

        let raw_floats: *mut f32 = unsafe {
            MLEM_TOF_Reco(
                niterations,
                TOF, TOF_resolution,
                FOV_XY, FOV_Z, NXY, NZ,
                ncoinc,
                LOR_X1, LOR_Y1, LOR_Z1, LOR_T1,
                LOR_X2, LOR_Y2, LOR_Z2, LOR_T2,
                SENS,
                outfile_prefix,
                out_niter,
            )
        };
        let image = unsafe {
            std::slice::from_raw_parts(raw_floats, image_size)
            .into_iter()
            .collect::<Vec<_>>()
        };
        assert_eq!(image.len(), image_size);
    }
}
