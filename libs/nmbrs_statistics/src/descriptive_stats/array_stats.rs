use crate::AsSlice;
use nmbrs_algebra::{MidPoint, NumericField, NumericSemiGroup};

/* TODOs:
- splt into descriptive and inferential stats and ordered and unordered stats
- math formulas
- combined stats (returning mean and std) for perf reasons
- furhter metrics */

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

/// Calculates the [empirical percentile](https://en.wikipedia.org/wiki/Percentile) of the _sorted_ samples.
/// The samples are assumed to be sorted in ascending order and level is assumed to be in the range `[0, 1]`.
pub fn percentile<T>(sorted_xs: &[T], level: f64) -> Option<T>
where
    T: NumericField + MidPoint + Copy,
{
    if !(0.0..=1.0).contains(&level) {
        return None;
    }
    if sorted_xs.is_empty() {
        return None;
    }

    let n = sorted_xs.len();

    // NOTE: have to add `-1` below due to (mathematical) idx start of 1 (rather than 0)
    let candidate_idx: f64 = n as f64 * level;
    let floored: usize = candidate_idx.floor() as usize;

    // case candidate is an integer
    if candidate_idx == floored as f64 {
        let idx_bottom = (floored - 1).max(0);
        let idx_top = floored.min(n);
        return Some(sorted_xs[idx_bottom].mid_point(sorted_xs[idx_top]));
    }
    let idx = ((candidate_idx + 1.0).floor().min(n as f64) as usize - 1).max(0);
    Some(sorted_xs[idx])
}

pub trait Percentile<T> {
    fn percentile(&self, level: f64) -> Option<T>;

    fn median(&self) -> Option<T> {
        self.percentile(0.5)
    }

    fn p75(&self) -> Option<T> {
        self.percentile(0.75)
    }

    fn p25(&self) -> Option<T> {
        self.percentile(0.25)
    }
}

impl<'a, T, S> Percentile<T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericField + MidPoint + Copy,
{
    fn percentile(&self, level: f64) -> Option<T> {
        percentile(self.as_slice(), level)
    }
}

#[cfg(test)]
mod test {
    use crate::descriptive_stats::array_stats::VarianceBias;

    use super::{Covariance, Dot, Mean, Percentile, Sum, Variance};

    #[test]
    fn sum() {
        assert_eq!(super::sum(&Vec::with_capacity(0)) as Option<f64>, None);
        assert_eq!(super::sum(&[]) as Option<f64>, None);

        let xs = vec![1.0, 1.0, 2.0];

        assert_eq!(super::sum(&xs), xs.sum());
        assert_eq!(super::sum(&xs), Some(4.0));

        let xs = vec![1.0, 2.0, 3.5];
        assert_eq!(super::sum(&xs), Some(6.5));
    }

    #[test]
    fn mean() {
        assert_eq!(super::mean(&Vec::with_capacity(0)) as Option<f64>, None);

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
            super::variance(&xs, Some(VarianceBias::Population)),
            Some(2.0)
        );

        let _xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // assert_eq!(super::variance(&xs), Some(2.));
        // assert_eq!(super::variance(&xs), super::covariance(&xs, &xs)); // TODO: rescale
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

    #[test]
    fn percentile() {
        let mut samples = vec![82., 91., 12., 92., 63., 9., 28., 55., 96., 97.];
        samples.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let median = super::percentile(&samples, 0.5);
        assert_eq!(median, Some(72.5));
        assert_eq!(super::percentile(&samples, 0.5), samples.percentile(0.5));
        assert_eq!(super::percentile(&samples, 0.5), samples.median());

        let quartile_fst = super::percentile(&samples, 0.25);
        assert_eq!(quartile_fst, Some(28.0));
        assert_eq!(super::percentile(&samples, 0.25), samples.p25());

        let quartile_trd = super::percentile(&samples, 0.75);
        assert_eq!(quartile_trd, Some(92.0));
        assert_eq!(super::percentile(&samples, 0.75), samples.p75());
    }

    use approx::assert_abs_diff_eq;
    const EPSILON: f64 = 1e-15;

    #[test]
    fn mean1() {
        let data = [
            5.3766713954610001e-01,
            1.8338850145950865e+00,
            -2.2588468610036481e+00,
            8.6217332036812055e-01,
            3.1876523985898081e-01,
            -1.3076882963052734e+00,
            -4.3359202230568356e-01,
            3.4262446653864992e-01,
            3.5783969397257605e+00,
            2.7694370298848772e+00,
        ];
        assert_abs_diff_eq!(
            super::mean(&data).unwrap(),
            6.2428219709029698e-01,
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
            5.3766713954610001e-01,
            1.8338850145950865e+00,
            -2.2588468610036481e+00,
            8.6217332036812055e-01,
            3.1876523985898081e-01,
            -1.3076882963052734e+00,
            -4.3359202230568356e-01,
            3.4262446653864992e-01,
            3.5783969397257605e+00,
            2.7694370298848772e+00,
        ];
        assert_abs_diff_eq!(
            super::variance(xs, Some(VarianceBias::Sample)).unwrap(),
            3.1324921339484746e+00,
            epsilon = EPSILON
        );
    }

    // #[test]
    // fn variance_0() {
    //     assert!(super::population_variance::<f64>(&[]).is_none());
    // }

    // #[test]
    // fn variance_1() {
    //     assert_eq!(super::population_variance::<f64>(&[1.0]), Some(0.0));
    // }
}
