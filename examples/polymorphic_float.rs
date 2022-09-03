/// This example demonstrates how to write code which is generic over floats
/// with different precisions. It depends on the `Float` and `FromPrimitive`
/// traits from the `num-traits` crate.
///
/// As a secondary concern, it compares execution speeds of f32 vs f64
/// operations.
///
/// It is a generalization of the `azip` example, but the emphasis here is on
/// float-polymorphism rather than serial vs. parallel iteration.

// This is the trait which allows us to abstract over different float types.
use num_traits::Float;

// This trait is needed in order to convert from primitive float types to
// generic ones. In this example we need it in order to specify specific float
// values by hand.
//
// NOTE: If such values get passed to some other polymorphic function, they
// taint it with this trait too.
use num_traits::FromPrimitive;

// A function which is polymorphic over float types. Note the use of helper
// function `f` (defined just below) to create literal values.
fn calculate<N: Float + FromPrimitive>(b: &N, c: &N) -> N {
    (b.sin().powf(f(2.1)) + c.cos().powf(f(1.9))).sqrt()
}

// Helper function to mitigate the boilerplate needed to convert from concrete
// float types to the generic one.
fn f<N: Float + FromPrimitive>(x: f32) -> N { N::from_f32(x).unwrap() }

// Another float-polymorphic function. This one processes the values on multiple
// threads (via par_azip), so we need to explicitly constrain our abstract type
// to `Send + Sync`.
fn timings<N: Float + FromPrimitive + Send + Sync>() -> Timings {

    type A<N> = Array3<N>;

    let mut t = Instant::now();

    let n = Cli::from_args();
    let mut a = A::<N>::zeros((n.x, n.y, n.z));
    let     b = A::from_elem(a.dim(), f(0.1));
    let     c = A::from_elem(a.dim(), f(0.2));
    let alloc = report_time(&mut t, "allocated");

    par_azip!((a in &mut a, &b in &b, &c in &c) *a = calculate(&b, &c));
    let par = report_time(&mut t, "parallel");

    azip!((a in &mut a, &b in &b, &c in &c) *a = calculate(&b, &c));
    let ser = report_time(&mut t, "serial");

    println!("Speedup factor {:?}", ratio(ser, par));

    Timings{ alloc, par, ser }
}

// Client code that uses the `timings` function above, with two different
// concrete float types.
fn main() {
    println!("--------- f32 ----------");
    let tf32 = timings::<f32>();
    println!("--------- f64 ----------");
    let tf64 = timings::<f64>();

    println!("----- f64 timings / f32 timings -----");
    println!("Allocation: {}", ratio(tf64.alloc, tf32.alloc));
    println!("Parallel  : {}", ratio(tf64.par  , tf32.par  ));
    println!("Serial    : {}", ratio(tf64.ser  , tf32.ser  ));
}

// ================================================================================
// ================= Details not directly relevant to the example =================
// ================================================================================
use ndarray::{Array3, azip, parallel::par_azip};
// -------------------- Utility for reporting timings --------------------
use std::time::{Duration, Instant};

fn report_time(t: &mut Instant, message: &str) -> Duration {
    let elapsed = t.elapsed();
    println!("{}: {} ms", message, elapsed.as_millis());
    *t = Instant::now();
    elapsed
}

struct Timings {
    alloc: Duration,
    par: Duration,
    ser: Duration,
}

fn ratio(n: Duration, d: Duration) -> f32 {
    n.as_secs_f32() / d.as_secs_f32()
}
// ----------------------------------- CLI -----------------------------------
use clap::Parser;

#[derive(Parser, Debug)]
#[structopt(setting = clap::AppSettings::ColoredHelp)]
#[structopt(name = "azip", about = "Compare parallel and sequential zipping in ndarray")]
pub struct Cli {

    /// Size of array in fastest-changing index dimension
    #[structopt(default_value = "256")]
    x: usize,

    /// Size of array in      middle            dimension
    #[structopt(default_value = "256")]
    y: usize,

    /// Size of array in slowest-changing index dimension
    #[structopt(default_value = "256")]
    z: usize,

}
