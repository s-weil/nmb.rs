mod samples;

use algebra::{NumericField, NumericSemiGroup};
use samples::AsSlice;

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
pub fn variance<T>(xs: &[T]) -> Option<T>
where
    T: NumericField + From<i8> + Copy,
{
    let len = xs.len() as i8;
    let mean = mean(xs)?;

    let mse = xs.iter().fold(T::zero(), |err, x| {
        let x_err = -*x + mean;
        err + x_err * x_err
    });

    Some(mse / len.into())
}

pub trait Variance<T> {
    fn variance(&self) -> Option<T>;
}

impl<'a, T, S> Variance<T> for S
where
    S: AsSlice<T>,
    T: 'a + NumericField + From<i8> + Copy,
{
    fn variance(&self) -> Option<T> {
        variance(self.as_slice())
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

#[cfg(test)]
mod test {
    use crate::{Covariance, Dot, Mean, Sum, Variance};

    #[test]
    fn sum() {
        assert_eq!(super::sum(&Vec::with_capacity(0)) as Option<f64>, None);
        assert_eq!(super::sum(&vec![]) as Option<f64>, None);

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
    fn variance() {
        let xs: Vec<f32> = Vec::with_capacity(0);
        assert_eq!(super::variance(&xs), None);

        let xs = vec![2.0, 2.0, 2.0, 2.0, 2.0];
        assert_eq!(super::variance(&xs), Some(0.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(super::variance(&xs), Some(2.0));
        assert_eq!(super::variance(&xs), xs.variance());

        let xs: &[f32] = &[1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(super::variance(&xs), Some(2.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
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

        // let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        // assert_eq!(super::covariance(&xs, &ys), Some(154.0));
    }
}
