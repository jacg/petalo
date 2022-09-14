use cpython::{Python, PyDict, PyResult};

fn main() {
    let gil = Python::acquire_gil();
    hello(gil.python()).unwrap();
}

fn hello(py: Python) -> PyResult<()> {
    let sys = py.import("sys")?;
    let version = get_python_version(py)?;

    let locals = PyDict::new(py);
    locals.set_item(py, "os", py.import("os")?)?;
    let user: String = py.eval("os.getenv('USER') or os.getenv('USERNAME')", None, Some(&locals))?.extract(py)?;

    println!("Hello {}, I'm Python {}", user, version);
    Ok(())
}

fn get_python_version(py: Python) -> PyResult<String> {
    let sys = py.import("sys")?;
    sys.get(py, "version")?.extract(py)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn python_version() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        assert_eq!(get_python_version(py)?, "3.10.4 (main, Mar 23 2022, 20:25:24) [GCC 11.3.0]");
        Ok(())
    }
}
