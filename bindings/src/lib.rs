use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use petalo::{mlem::Image, weights::VoxelBox};

#[pyfunction]
/// Calculate CRC for a 60x60x60 voxel image
fn crcs(data: Vec<f32>) -> Vec<f32> {

    // Regions of interest for CRC
    fn polar(r: f32, phi: f32) -> (f32, f32) { (r * phi.cos(), r * phi.sin()) }

    use petalo::fom::ROI;
    let step = std::f32::consts::PI / 6.0;
    let roi_from_centre = 50.0;
    let (hot, cold, bg_activity, bg_radius) = (4.0, 0.0, 1.0, 4.0);
    let rois = vec![
        (ROI::CylinderZ(polar(roi_from_centre,  2.0*step),  4.0),  hot),
        (ROI::CylinderZ(polar(roi_from_centre,  4.0*step),  6.5),  hot),
        (ROI::CylinderZ(polar(roi_from_centre,  6.0*step),  8.5),  hot),
        (ROI::CylinderZ(polar(roi_from_centre,  8.0*step), 11.0),  hot),
        (ROI::CylinderZ(polar(roi_from_centre, 10.0*step), 14.0), cold),
        (ROI::CylinderZ(polar(roi_from_centre, 12.0*step), 18.5), cold),
    ];

    let bg_rois = vec![
        ROI::CylinderZ(polar(roi_from_centre,  1.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  3.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  5.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  7.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre,  9.0*step), bg_radius),
        ROI::CylinderZ(polar(roi_from_centre, 11.0*step), bg_radius),
    ];

    let vbox = VoxelBox::new((180.0, 180.0, 180.0), (60, 60, 60));
    let image = Image::new(vbox, data);
    let crcs = image.foms(&rois, &bg_rois, bg_activity, true).crcs;
    crcs
}

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
    m.add_function(wrap_pyfunction!(crcs, m)?)?;
    m.add_function(wrap_pyfunction!(rust_enum_parameter, m)?)?;

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

    Ok(())
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
