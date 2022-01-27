use ndhistogram::axis::{Axis, BinInterval};
use num_traits::{Float, Num, NumCast, NumOps};
use serde::{Deserialize, Serialize};

/// A wrap-around axis with equal-sized bins.
///
/// An axis with `N` equally-spaced, equal-sized bins, in `[low, high)`.
/// Entries outside this interval get wrapped around.
/// There are no overflow bins so this axis has exactly `N` bins.
///
/// # Examples
/// 1D histogram representing TODO
/// ```
/// // use ndhistogram::{ndhistogram, Histogram};
/// // use ndhistogram::axis::{Axis, BinInterval};
/// // use petalo::sinogram::Cyclic;
/// // let mut hist = ndhistogram!(Cyclic::new(4, 0.0, 360.0));
/// // let axis = &hist.axes().as_tuple().0;
/// // hist.
/// ```
#[derive(Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub struct Cyclic<T = f64> {
    nbins: usize,
    low: T,
    high: T,
    step: T,
}

impl<T> Cyclic<T>
where
    T: PartialOrd + Num + NumCast + NumOps + Copy,
{
    /// Create a wrap-around axis with `nbins` uniformly-spaced bins in the range `[low, high)`.
    ///
    /// Only implemented for [Float]. TODO Use [Cyclic::with_step_size] for integers.
    ///
    /// # Panics
    /// Panics if `nbins == 0` or `low == high`.
    pub fn new(nbins: usize, low: T, high: T) -> Self
    where
        T: Float
    {
        if nbins == 0  { panic!("Need more than zero bins on axis") }
        if low == high { panic!("Axis range must be non-zero") }
        let (low, high) = if low < high { (low, high) } else { (high, low) };
        let step = (high - low) / T::from(nbins).unwrap();
        Self { nbins, low, high, step }
    }

    /// Create a wrap-around axis with `nbins` uniformly-spaced bins in the range `[low, low+num*step)`.
    /// # Panics
    /// Panics if `nbins == 0` or `step <= 0`.
    pub fn with_step_size(nbins: usize, low: T, step: T) -> Self {
        let high = T::from(nbins).expect("Failed to convert nbins to coordinate type") * step + low;
        if nbins == 0        { panic!("Need more than zero bins on axis") }
        if step <= T::zero() { panic!("Step size must be strictly positive") }
        Self { nbins, low, high, step }
    }
}

impl<T> Cyclic<T> {
    /// Low edge of axis (excluding wrap-around) // TODO or should this be - infinity?
    pub fn low(&self) -> &T { &self.low }
    /// High edge of axis (excluding wrap-around) // TODO or should this be + infinity?
    pub fn high(&self) -> &T { &self.high }
}


// TODO integers?
impl<T: PartialOrd + NumCast + NumOps + Copy> Axis for Cyclic<T> {
    type Coordinate = T;
    type BinInterval = BinInterval<T>;

    #[inline]
    fn index(&self, coordinate: &Self::Coordinate) -> Option<usize> {
        let steps = (*coordinate - self.low) / self.step;

        let c = coordinate.to_f32().unwrap();
        let l = self.low.to_f32().unwrap();
        let s = self.step.to_f32().unwrap();
        let st = steps.to_f32().unwrap();
        println!("({c} - {l}) / {s} = {st}");

        Some((steps.to_usize().expect("TODO")) % self.nbins)
    }

    fn num_bins(&self) -> usize { self.nbins }

    fn bin(&self, index: usize) -> Option<<Self as Axis>::BinInterval> {
        todo!()
    }
}

#[cfg(test)]
mod test_index {
    use ndhistogram::axis::Uniform;
    use rstest::rstest;

    use super::*;

    #[rstest( bin_no, expected_interval,
              case(0, Some(BinInterval::new(0.00, 0.25))),
              case(1, Some(BinInterval::new(0.25, 0.50))),
              case(2, Some(BinInterval::new(0.50, 0.75))),
              case(3, Some(BinInterval::new(0.75, 1.00))),
    )]
    fn bin(bin_no: usize, expected_interval: Option<BinInterval<f32>>) {
        let axis = Cyclic::new(4, 0.0, 1.0);
        assert_eq!(axis.bin(bin_no), expected_interval);
    }

    #[rstest( bin_no, expected_interval,
              case(0, Some(BinInterval::underflow(0.0))),
              case(1, Some(BinInterval::new(0.00, 0.25))),
              case(0, Some(BinInterval::underflow(0.0))),
              case(1, Some(BinInterval::new(0.00, 0.25))),
              case(2, Some(BinInterval::new(0.25, 0.50))),
              case(3, Some(BinInterval::new(0.50, 0.75))),
              case(4, Some(BinInterval::new(0.75, 1.00))),
              case(5, Some(BinInterval::overflow(1.0))),
    )]
    fn bin_uniform(bin_no: usize, expected_interval: Option<BinInterval<f32>>) {
        let axis = Uniform::new(4, 0.0, 1.0);
        assert_eq!(axis.bin(bin_no), expected_interval);
    }

    #[rstest(/**/ coordinate,      expected_index,
             case(  0.0 , Some(0)),
             case(  0.09, Some(0)),
             case(  0.1 , Some(1)),
             case(  0.19, Some(1)),
             case(  0.2 , Some(2)),
             case( 10.0 , Some(0)),
             case( 20.33, Some(3)),
             case( 50.99, Some(9)),
             case( -0.1 , Some(9)),
             case( -0.19, Some(9)),
             case( -0.2 , Some(8)),
             case( -0.9 , Some(1)),
             case( -0.95, Some(0)),
             case(-10.0 , Some(0)),
             case(-10.05, Some(9)),
             case(-10.1 , Some(8)),
    )]
    fn index(coordinate: f32, expected_index: Option<usize>) {
        let axis = Cyclic::new(10, 0.0, 1.0);
        assert_eq!(axis.index(&coordinate), expected_index);
    }

    #[test]
    fn indices() {
        let n = 7;
        let axis = Cyclic::new(n, 23.4, 97.3);
        let indices = axis.indices().collect::<Vec<_>>();
        assert_eq!(indices, (0..n).collect::<Vec<_>>());
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn wrap() {
        use super::*;
        use ndhistogram::{ndhistogram};
        let mut hist = ndhistogram!(Cyclic::new(4, 0.0, 360.0));


    }
}
