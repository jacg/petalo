/// This example shows how to perform an elementwise operation on three ndarrays
/// by zipping over them, both sequentially and in parallel.
///
/// It reports and compares the times taken by both approaches.
///
/// The parallel version requires the `rayon` feature of `ndarray` to be
/// enabled.
///
/// To run this example:
///
///    cargo run --example azip --release
///
/// If you want to control the size of the array:
///
///    cargo run --example azip --release 100         # size: (100 * 256 * 256)
///    cargo run --example azip --release 100 200     # size: (100 * 200 * 256)
///    cargo run --example azip --release 100 200 300 # size: (100 * 200 * 300)

// The serial version of the macro which abstracts zipping over multiple arrays
// of the same shape.
use ndarray::azip;
// The parallel version of azip. Requires the "rayon" feature of ndarray to be
// enabled.
use ndarray::parallel::par_azip;

// Compare the times taken
fn main() {

    let mut t = Instant::now();
    let n = Cli::from_args();

    // Allocation of arrays
    let mut a = A::zeros((n.x, n.y, n.z));
    let     b = A::from_elem(a.dim(), 1.);
    let     c = A::from_elem(a.dim(), 2.);
    report_time(&mut t, "allocation");

    // Parallel calculation
    par_azip!((a in &mut a, &b in &b, &c in &c) *a = calculate(&b, &c));
    let par_time = report_time(&mut t, "parallel");

    // Sequential calculation. Only difference: azip! -> par_azip!
    azip!((a in &mut a, &b in &b, &c in &c) *a = calculate(&b, &c));
    let ser_time = report_time(&mut t, "serial");

    // Report speedup
    println!("Speedup factor: (serial time) / (parallel time) {:?}",
             ser_time.as_secs_f32() / par_time.as_secs_f32());
}

// An arbitrary calculation that burns some CPU cycles. It acts on a single set
// of homologous values in the array, at a time.
fn calculate(b: &N, c: &N) -> N {
    (b.sin().powf(2.1) + c.cos().powf(1.9)).sqrt()
}


type N = f32;

use ndarray::Array3;
type A = Array3<N>;

// -------------------- Utility for reporting timings --------------------
use std::time::Instant;

fn report_time(t: &mut Instant, message: &str) -> std::time::Duration {
    let elapsed = t.elapsed();
    println!("{}: {} ms", message, elapsed.as_millis());
    *t = Instant::now();
    elapsed
}

// ----------------------------------- CLI -----------------------------------
use clap::Parser;

#[derive(Parser, Debug)]
#[clap(setting = clap::AppSettings::ColoredHelp)]
#[clap(name = "azip", about = "Compare parallel and sequential zipping in ndarray")]
pub struct Cli {

    /// Size of array in fastest-changing index dimension
    #[clap(default_value = "256")]
    x: usize,

    /// Size of array in      middle            dimension
    #[clap(default_value = "256")]
    y: usize,

    /// Size of array in slowest-changing index dimension
    #[clap(default_value = "256")]
    z: usize,

}
