pub mod axis;

use ndhistogram::{ndhistogram, VecHistogram, axis::Uniform, AxesTuple, HistND};

use std::f32::consts::{PI,TAU, FRAC_PI_2 as PI_2};

// type HistND<A, V = f64> = VecHistogram<AxesTuple<A>, V>;

type Hmm = HistND<(Uniform<f32>,
                   Uniform<f32>,
                   Uniform<f32>,
                   Uniform<f32>,),
                  usize>;


pub struct Sinogram {
    histogram: Hmm,
}

impl Sinogram {
    fn new(r: f32, l: f32) -> Self {
        Self {
            histogram: ndhistogram!(
                Uniform::new(10, -l/2.0, l/2.0), // z
                Uniform::new(10,    0.0, r    ), // r
                // TODO remove overflow bins for angles, create own wrap around axis
                Uniform::new(10,    0.0, TAU  ), // phi  : around z, anticlockwise from x
                Uniform::new(10,  -PI_2, PI_2 ); // theta: around r, anticlockwise from perpendicular to z
                usize
            )
        }
    }
}

//type Histywisty = VecHistogram<>;

#[cfg(test)]
mod test {

    #[test]
    fn foo() {
        assert!(2 == 8 / 4);
    }

    #[test]
    fn bar() {
        assert!(2 == 8 / 4);
    }

}
