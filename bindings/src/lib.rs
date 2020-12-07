use pyo3::prelude::*;
use pyo3::wrap_pyfunction;


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
