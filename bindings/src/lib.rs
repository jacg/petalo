use pyo3::prelude::*;
use pyo3::wrap_pyfunction;


#[pyfunction]
/// The naive, recursive fibonacci implementation
fn fib(n: usize) -> usize {
    if n < 2 { 1 }
    else     { fib(n-1) + fib(n-2) }
}

#[pymodule]
fn fulano(_python_gil: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fib, m)?)?;
    Ok(())
}
