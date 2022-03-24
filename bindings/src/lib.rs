use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use petalo::{mlem::Image, weights::FOV, fom, types::{Length as L, Intensity}};

#[pyfunction]
#[text_signature = "(n, /)"]
/// The naive, recursive fibonacci implementation
fn fib(n: usize) -> usize {
    if n < 2 { 1 }
    else     { fib(n-1) + fib(n-2) }
}

#[pymodule]
/// Module docstring works too!
fn fulano(_py_gil: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fib, m)?)?;
    m.add_function(wrap_pyfunction!(rust_enum_parameter, m)?)?;
    m.add_function(wrap_pyfunction!(roi, m)?)?;
    m.add_function(wrap_pyfunction!(fom_config, m)?)?;

    #[pyfn(m, "fab")]
    #[text_signature = "(n, /)"]
    /// The iterative fibonacci implementation
    fn burp(_py_gil: Python, mut n: usize) -> usize {
        let (mut p, mut c) = (0,1);
        while n > 0 {
            let old_c = c;
            c = c + p;
            p = old_c;
            n -= 1;
        }
        c
    }

    m.add_class::<Lift>()?;
    m.add_class::<FomConfig>()?;

    Ok(())
}

#[pyclass]
struct FomConfig {
    cfg: fom::FomConfig,
    fov: FOV,
}

#[pymethods]
impl FomConfig {

    #[new]
    fn new(rois: Vec<(ROI, Intensity)>, bg_rois: Vec<ROI>, bg: Intensity, voxels: (usize, usize, usize), size: (L,L,L)) -> Self {
        let rois: Vec<(petalo::fom::ROI, Intensity)> = rois.into_iter()
            .map(|(r,i)| (pyroi_to_fomroi(r), i))
            .collect();
        let background_rois = bg_rois.into_iter().map(pyroi_to_fomroi).collect();

        let cfg = fom::FomConfig{ rois, background_rois, background_activity: bg};
        FomConfig{ cfg, fov: FOV::new(size, voxels)}
    }

    /// Calculate CRC for a 60x60x60 voxel image
    fn crcs(&self, data: Vec<Intensity>) -> Vec<Intensity> {
        let image = Image::new(self.fov, data);
        let crcs = image.foms(&self.cfg, true).crcs;
        crcs
    }

}


#[pyclass]
#[text_signature = "(initial_height)"]
/// It's a Lift: it goes up and down
struct Lift {
    #[pyo3(get)]
    height: i32
}

#[pymethods]
impl Lift {

    #[new] // Signature goes on the struct
    fn new(initial_height: i32) -> Self { Self { height: initial_height }}

    fn up  (&mut self, n: usize) { self.height += n as i32 }
    fn down(&mut self, n: usize) { self.height -= n as i32 }

}

#[pyfunction]
fn roi(roi: ROI) -> String {
    use ROI::*;
    match roi {
        Sphere{x, y, z, r} => format!("S {} {} {} {}", x, y, z, r),
        CylinderX{y, z, r} => format!("X {} {} {}", y, z, r),
        CylinderY{x, z, r} => format!("Y {} {} {}", x, z, r),
        CylinderZ{x, y, r} => format!("Z {} {} {}", x, y, r),
    }
}


fn pyroi_to_fomroi(pyroi: ROI) -> petalo::fom::ROI {
    use              ROI as lr;
    use petalo::fom::ROI as fr;
    match pyroi {
        lr::Sphere {x,y,z,r} => fr::Sphere ((x,y,z), r),
        lr::CylinderX{y,z,r} => fr::CylinderX((y,z), r),
        lr::CylinderY{x,z,r} => fr::CylinderX((x,z), r),
        lr::CylinderZ{x,y,r} => fr::CylinderZ((x,y), r),
    }
}

#[pyfunction]
fn fom_config(rois: Vec<(ROI, Intensity)>, bg_rois: Vec<ROI>, bg: Intensity) -> String /*FomConfig*/ {
    let rois: Vec<(petalo::fom::ROI, Intensity)> = rois.into_iter()
        .map(|(r,i)| (pyroi_to_fomroi(r), i))
        .collect();
    let background_rois = bg_rois.into_iter().map(pyroi_to_fomroi).collect();

    let config = fom::FomConfig{ rois, background_rois, background_activity: bg};
    format!("{:?}", config)
}

#[derive(FromPyObject)]
enum ROI {
    Sphere{ x: L, y: L, z: L, r: L },
    CylinderZ{ x: L, y: L, r: L },
    CylinderY{ x: L, z: L, r: L },
    CylinderX{ y: L, z: L, r: L },
}


#[pyfunction]
/// Testing Rust enum conversion
fn rust_enum_parameter(e: RustyEnum) -> String {
    use RustyEnum::*;
    match e {
        Int(n)                   => format!("Int({})", n),
        String(s)                => format!("String(\"{}\")", s),
        IntTuple(a,b)            => format!("IntTuple({}, {})", a, b),
        StringIntTuple(a,b)      => format!("StringTuple(\"{}\", {})", a, b),
        Coordinates3d {x, y, z}  => format!("Coordinates3d({}, {}, {})", x,y,z),
        Coordinates2d {a:x, b:y} => format!("Coordinates2d({}, {})"    , x,y),
        //CatchAll(pyany)          => format!("CatchAll: {:?}", pyany),
    }
}

#[derive(FromPyObject)]
enum RustyEnum {
    Int(usize), // input is a positive int
    String(String), // input is a string
    IntTuple(usize, usize), // input is a 2-tuple with positive ints
    StringIntTuple(String, usize), // input is a 2-tuple with String and int
    Coordinates3d { // needs to be in front of 2d
        x: usize,
        y: usize,
        z: usize,
    },
    Coordinates2d { // only gets checked if the input did not have `z`
        #[pyo3(attribute("x"))]
        a: usize,
        #[pyo3(attribute("y"))]
        b: usize,
    },
    //#[pyo3(transparent)]
    //CatchAll(&'a PyAny), // This extraction never fails
}
