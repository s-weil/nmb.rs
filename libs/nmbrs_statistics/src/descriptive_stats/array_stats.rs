use crate::AsSlice;
use nmbrs_algebra::{NumericField, NumericSemiGroup};

/*
Array statistics provides routines optimized for single-dimensional arrays.
 */

// TODO: don't need a ring here
pub fn sum<T>(xs: &[T]) -> Option<T>
where
    T: NumericSemiGroup + Copy,
{
    if xs.is_empty() {
        return None;
    }

    let sum = xs.iter().fold(T::zero(), |acc, x| acc + *x);
    Some(sum)
}

pub trait Sum<T> {
    fn sum(&self) -> Option<T>;
}

// impl<'a, T, S> Sum<T> for S
// where
//     S: AsRef<&'a [T]> + ?Sized,
//     T: 'a + NumericRing + Copy,
// {
//     fn sum(&self) -> Option<T> {
//         sum2(self.as_ref())
//     }
// }

impl<'a, T, S> Sum<T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericSemiGroup + Copy,
{
    fn sum(&self) -> Option<T> {
        sum(self.as_slice())
    }
}

/// The arithmetic mean or average of the provided samples.
/// In statistics, the sample mean is a measure of the central tendency and estimates the expected value of the distribution.
/// The mean is affected by outliers, so if you need a more robust estimate consider to use the Median instead.
/// $ \Sigma x_i $
/// Sources:
/// - [MathDotNet](https://numerics.mathdotnet.com/DescriptiveStatistics)
/// - [Wikipedia](https://en.wikipedia.org/wiki/Mean)
/// - [Wolfram MathWorld](http://mathworld.wolfram.com/SampleMean.html)
pub fn mean<T>(xs: &[T]) -> Option<T>
where
    T: NumericField + From<i8> + Copy,
{
    let len = xs.len() as i8;
    let sum: T = sum(xs)?;

    Some(sum / len.into())
}

pub trait Mean<T> {
    fn mean(&self) -> Option<T>;
}

impl<'a, T, S> Mean<T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericField + From<i8> + Copy,
{
    fn mean(&self) -> Option<T> {
        mean(self.as_slice())
    }
}

/// The (biased) [sample variance](https://en.wikipedia.org/wiki/Variance#Sample_variance).
///
/// NOTE: The variance is covered by the `Covariance` but provided as a more performant function.
/// Sources:
/// * [MathDotNet](https://numerics.mathdotnet.com/DescriptiveStatistics)
/// * [Wikipedia](https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Two-pass_algorithm)
pub fn variance<T>(xs: &[T], ty: Option<VarianceBias>) -> Option<T>
where
    T: NumericField + From<i8> + Copy,
{
    let len = xs.len() as i8;
    let mean = mean(xs)?;

    let mse = xs.iter().fold(T::zero(), |err, x| {
        let x_err = -*x + mean;
        err + x_err * x_err
    });

    let scale = ty.unwrap_or_default().scale(len);
    Some(mse / scale.into())
}

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarianceBias {
    /// Biased estimator of the population variance.
    Population,
    /// Unbiased estimator for the sample variance.
    #[default]
    Sample,
}

impl VarianceBias {
    fn scale(&self, n: i8) -> i8 {
        match self {
            VarianceBias::Population => n,
            VarianceBias::Sample => n - 1,
        }
    }
}

pub trait Variance<T> {
    /// TODO: add docs as above
    fn sample_variance(&self) -> Option<T>;
    fn population_variance(&self) -> Option<T>;
}

impl<'a, T, S> Variance<T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericField + From<i8> + Copy,
{
    fn population_variance(&self) -> Option<T> {
        variance(self.as_slice(), Some(VarianceBias::Population))
    }
    fn sample_variance(&self) -> Option<T> {
        variance(self.as_slice(), Some(VarianceBias::Sample))
    }
}

pub fn dot<T>(xs: &[T], ys: &[T]) -> Option<T>
where
    T: NumericField + From<i8> + Copy,
{
    if xs.is_empty() || xs.len() != ys.len() {
        return None;
    }

    let dot = xs
        .iter()
        .zip(ys.iter())
        .fold(T::zero(), |acc, (x, y)| acc + *x * *y);
    Some(dot)
}

pub trait Dot<S, T> {
    fn dot(&self, ys: S) -> Option<T>;
}

impl<'a, T, S> Dot<S, T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericField + From<i8> + Copy,
{
    fn dot(&self, ys: S) -> Option<T> {
        dot(self.as_slice(), ys.as_slice())
    }
}
// TODO: create a macro for dot

/// https://en.wikipedia.org/wiki/Sample_mean_and_covariance
pub fn covariance<T>(xs: &[T], ys: &[T]) -> Option<T>
where
    T: NumericField + From<i8> + Copy,
{
    if xs.len() != ys.len() || xs.len() <= 1 {
        return None;
    }

    let x_mean = mean(xs)?;
    let y_mean = mean(ys)?;

    // TODO: bench it. maybe it's faster to use a single loop
    let x_err: Vec<T> = xs.iter().map(|x| *x - x_mean).collect();
    let y_err: Vec<T> = ys.iter().map(|y| *y - y_mean).collect();

    let dot = dot(&x_err, &y_err)?;
    let len = xs.len() as i8;
    Some(dot / (T::from(len) - T::one()))
}

pub trait Covariance<S, T> {
    fn covariance(&self, ys: S) -> Option<T>;
}

impl<'a, T, S> Covariance<S, T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericField + From<i8> + Copy,
{
    fn covariance(&self, ys: S) -> Option<T> {
        covariance(self.as_slice(), ys.as_slice())
    }
}

// TODO: add skewness
// https://en.wikipedia.org/wiki/Skewness#Sample_skewness

#[cfg(test)]
mod test {
    use super::{Covariance, Dot, Mean, Sum, Variance, VarianceBias};

    #[test]
    fn sum() {
        assert_eq!(super::sum(&Vec::with_capacity(0)) as Option<f64>, None);
        assert_eq!(super::sum::<i32>(&[]), None);

        let xs = [1, 2, 3, 4, 5, 6, 7, 8, 9];
        assert_eq!(super::sum::<i8>(&xs), Some(4 * 10 + 5));

        let xs = vec![1.0, 1.0, 2.0];
        assert_eq!(super::sum(&xs), xs.sum());
        assert_eq!(super::sum(&xs), Some(4.0));

        let xs = vec![1.0, 2.0, 3.5];
        assert_eq!(super::sum(&xs), Some(6.5));
    }

    #[test]
    fn mean() {
        assert_eq!(super::mean(&Vec::with_capacity(0)) as Option<f64>, None);
        assert!(super::mean::<i8>(&[]).is_none());

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(super::mean(&xs), xs.mean());
        assert_eq!(super::mean(&xs), Some(3.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(super::mean(&xs), xs.mean());
        assert_eq!(super::mean(&xs), Some(3.5));
    }

    #[test]
    fn population_variance() {
        let xs: Vec<f32> = Vec::with_capacity(0);
        assert_eq!(super::variance(&xs, Some(VarianceBias::Population)), None);

        let xs = vec![2.0, 2.0, 2.0, 2.0, 2.0];
        assert_eq!(
            super::variance(&xs, Some(VarianceBias::Population)),
            Some(0.0)
        );

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(
            super::variance(&xs, Some(VarianceBias::Population)),
            Some(2.0)
        );
        assert_eq!(
            super::variance(&xs, Some(VarianceBias::Population)),
            xs.population_variance()
        );

        let xs: &[f32] = &[1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(
            super::variance(xs, Some(VarianceBias::Population)),
            Some(2.0)
        );

        let _xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // assert_eq!(super::variance(&xs), Some(2.));
    }

    #[test]
    fn dot() {
        let xs = vec![1.0];
        let ys = vec![];
        assert_eq!(super::dot(&xs, &ys), None);

        let xs = vec![1.0];
        let ys = vec![4.0, 5.0];
        assert_eq!(super::dot(&xs, &ys), None);

        let xs = vec![1.0, 2.0, 3.0];
        let ys = vec![4.0, 5.0, 6.0];
        assert_eq!(super::dot(&xs, &ys), Some(32.0));
        assert_eq!(super::dot(&xs, &ys), xs.dot(ys));

        let xs = vec![1.0, 2.0, 3.0, 4.0];
        let ys = vec![4.0, 5.0, 6.0, 7.0];
        assert_eq!(super::dot(&xs, &ys), Some(60.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0];
        assert_eq!(super::dot(&xs, &ys), Some(100.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        assert_eq!(super::dot(&xs, &ys), Some(154.0));
    }

    #[test]
    fn covariance() {
        let xs = vec![1.0];
        let ys = vec![];
        assert_eq!(super::covariance(&xs, &ys), None);

        let xs = vec![1.0];
        let ys = vec![4.0, 5.0];
        assert_eq!(super::covariance(&xs, &ys), None);

        let xs = vec![1.0, 2.0, 3.0];
        let ys = vec![4.0, 5.0, 6.0];
        assert_eq!(super::covariance(&xs, &ys), Some(1.0));

        // let xs = vec![1.0, 2.0, 3.0, 4.0];
        // let ys = vec![4.0, 5.0, 6.0, 7.0];
        // assert_eq!(super::covariance(&xs, &ys), Some(60.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0];
        assert_eq!(super::covariance(&xs, &ys), Some(2.5));
        assert_eq!(super::covariance(&xs, &ys), xs.covariance(ys));

        assert_eq!(super::covariance(&xs, &xs), xs.sample_variance());

        // let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        // assert_eq!(super::covariance(&xs, &ys), Some(154.0));
    }

    use approx::assert_abs_diff_eq;
    const EPSILON: f64 = 1e-15;

    #[test]
    fn mean1() {
        let data = [
            5.376_671_395_461e-1,
            1.833_885_014_595_086_5,
            -2.258_846_861_003_648,
            8.621_733_203_681_206e-1,
            3.187_652_398_589_808e-1,
            -1.307_688_296_305_273_4,
            -4.335_920_223_056_835_6e-1,
            3.426_244_665_386_499e-1,
            3.578_396_939_725_760_5,
            2.769_437_029_884_877,
        ];
        assert_abs_diff_eq!(
            super::mean(&data).unwrap(),
            6.242_821_970_902_97e-1,
            epsilon = EPSILON
        );
    }

    #[test]
    fn mean_0() {
        assert_eq!(super::mean::<f64>(&[]), None);
    }

    #[test]
    fn sample_variance1() {
        let xs = &[
            5.376_671_395_461e-1,
            1.833_885_014_595_086_5,
            -2.258_846_861_003_648,
            8.621_733_203_681_206e-1,
            3.187_652_398_589_808e-1,
            -1.307_688_296_305_273_4,
            -4.335_920_223_056_835_6e-1,
            3.426_244_665_386_499e-1,
            3.578_396_939_725_760_5,
            2.769_437_029_884_877,
        ];
        assert_abs_diff_eq!(
            super::variance(xs, Some(VarianceBias::Sample)).unwrap(),
            3.132_492_133_948_474_6,
            epsilon = EPSILON
        );
    }
}
