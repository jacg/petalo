#![allow(nonstandard_style)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

use std::ffi::CString;
use ::std::os::raw::c_int;

pub fn cmlem(
    niterations: usize,
    TOF: bool,
    TOF_resolution: f32,
    FOV_XY: f32,
    FOV_Z : f32,
    NXY: usize,
    NZ : usize,
    ncoinc: usize,
    mut LOR_X1: Vec<f32>, mut LOR_Y1: Vec<f32>, mut LOR_Z1: Vec<f32>, mut LOR_T1: Vec<f32>,
    mut LOR_X2: Vec<f32>, mut LOR_Y2: Vec<f32>, mut LOR_Z2: Vec<f32>, mut LOR_T2: Vec<f32>,
    mut SENS: Vec<f32>,
    outfile_prefix: String,
    save_every_n: usize,
) {
    let outfile_prefix = CString::new(outfile_prefix).unwrap();
    unsafe {
        MLEM_TOF_Reco(
            niterations as c_int,
            TOF,
            TOF_resolution,
            FOV_XY,
            FOV_Z,
            NXY as c_int,
            NZ  as c_int,
            ncoinc as c_int,
            LOR_X1.as_mut_ptr(), LOR_Y1.as_mut_ptr(), LOR_Z1.as_mut_ptr(), LOR_T1.as_mut_ptr(),
            LOR_X2.as_mut_ptr(), LOR_Y2.as_mut_ptr(), LOR_Z2.as_mut_ptr(), LOR_T2.as_mut_ptr(),
            SENS.as_mut_ptr(),
            outfile_prefix.as_ptr(),
            save_every_n as c_int,
        );
    }
}

// TODO: write a test that uses the wrapper funcion. The only test we have so
// far uses the C version directly.

#[cfg(test)]
mod test {

    use super::*;
    use ::std::os::raw::c_char;

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

pub fn ale_burdel(dist: f32, deltaT: f32, TOF_resolution: f32) -> f32 {
    unsafe { ToFFunction(dist, deltaT, TOF_resolution) }
}

#[test]
fn test_burdel() {
    let burdel = ale_burdel(0.0, 0.0, 3.0);
    assert_eq!(burdel, 0.8871535);
}
