use ndhistogram::axis::{Axis, Uniform};
use ndhistogram::{ndhistogram, Histogram};

fn main() {
    let mut hist1d = ndhistogram!(Uniform::new(10, -5.0, 5.0));
    hist1d.fill(&-4.5);
    hist1d.fill(&-3.5);
    hist1d.fill(&-3.5);
    hist1d.fill_with(&-2.5, 2.1);
    hist1d.fill_with(&-1.5, 0.7);
    hist1d.fill_with(&-1.5, 1.7);
    hist1d.fill(&-100.0);
    //hist1d.fill_with_weighted(&-1.3, 2.3, 3);
    println!("{}", hist1d);

    let it = hist1d.value(&-2.9).unwrap();
    println!("{it} {it}");

    let nx = 20;
    let ny = 10;
    let mut hist2d = ndhistogram!(Uniform::new(nx, 0.0, 10.0),
                                  Uniform::new(ny, 0.0, 10.0));
    assert_eq!(hist2d.axes().num_bins(), (nx+2) * (ny+2));
    hist2d.fill(&(2.3,3.2));
    println!("{}", 8/2);

}
