use std::ops::{Add, Div, Mul, Neg};

pub trait Samples<'a, T: 'a>: std::iter::IntoIterator<Item = &'a T> {
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
    fn len(&self) -> usize;
}

// impl<'a, T> Samples<'a, T> for Vec<T> {
//     fn len(&self) -> usize {
//         (self as &Vec<T>).len()
//     }
// }

impl<'a, T> Samples<'a, T> for &'a Vec<T>
where
    T: Copy,
{
    fn len(&self) -> usize {
        (self as &'a Vec<T>).len()
    }
}

impl<'a, T> Samples<'a, T> for &'a [T] {
    fn len(&self) -> usize {
        (self as &[T]).len()
    }
}

/// Mimic some features of a ring, without ring axioms (Abelian group, Associativity, Distributivity)
pub trait NumericRing:
    Sized + AddIdentity + Add<Output = Self> + Mul<Output = Self> + Neg<Output = Self>
{
}

pub trait MulIdentity: Sized {
    fn one() -> Self;
}

pub trait NumericField: NumericRing + MulIdentity + Div<Output = Self> //Option<Self>>
{
}

pub trait AddIdentity: Sized {
    fn zero() -> Self;
}

#[macro_export]
macro_rules! impl_add_identity {
    ($impl_type:ty) => {
        impl AddIdentity for $impl_type {
            fn zero() -> Self {
                0 as $impl_type
            }
        }
    };
}

#[macro_export]
macro_rules! impl_mul_identity {
    ($impl_type:ty) => {
        impl MulIdentity for $impl_type {
            fn one() -> Self {
                1 as $impl_type
            }
        }
    };
}

// implement MulIdentity
impl_mul_identity! { i8 }
impl_mul_identity! { i16 }
impl_mul_identity! { i32 }
impl_mul_identity! { i64 }
impl_mul_identity! { f32 }
impl_mul_identity! { f64 }

// implement AddIdentity
impl_add_identity! { i8 }
impl_add_identity! { i16 }
impl_add_identity! { i32 }
impl_add_identity! { i64 }
impl_add_identity! { f32 }
impl_add_identity! { f64 }

impl NumericRing for f32 {}
impl NumericRing for f64 {}
impl NumericField for f32 {}
impl NumericField for f64 {}

// Todo: don't need a ring here
pub fn sum<'a, S, T: 'a>(xs: S) -> Option<T>
where
    S: Samples<'a, T>,
    T: NumericRing + Copy,
{
    if xs.is_empty() {
        return None;
    }

    let sum = xs.into_iter().fold(T::zero(), |acc, x| acc + *x);
    Some(sum)
}

pub trait Sum<T> {
    fn sum(self) -> Option<T>;
}

impl<'a, S, N: 'a> Sum<N> for S
where
    S: Samples<'a, N>,
    N: NumericRing + Copy,
{
    fn sum(self) -> Option<N> {
        sum(self)
    }
}

pub fn mean<'a, S, T: 'a>(xs: S) -> Option<T>
where
    S: Samples<'a, T>,
    T: NumericField + From<i8> + Copy,
{
    let len = xs.len() as i8;
    let sum: T = sum(xs)?;

    Some(sum / len.into())
}

pub trait Mean<T> {
    fn mean(self) -> Option<T>;
}

impl<'a, S, N: 'a> Mean<N> for S
where
    S: Samples<'a, N>,
    N: NumericField + From<i8> + Copy,
{
    fn mean(self) -> Option<N> {
        mean(self)
    }
}

/*



/// The (biased) sample variance: https://en.wikipedia.org/wiki/Variance#Sample_variance
/// NOTE: The variance is covered by the `Covariance` but provided as a more performant function.
pub fn variance(xs: &[f64]) -> Option<f64> {
    let mean = mean(xs)?;

    let mse = xs.iter().fold(0.0, |err, x| err + (x - mean).powi(2));
    Some(mse / xs.len() as f64)
}

pub fn dot(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.is_empty() || xs.len() != ys.len() {
        return None;
    }

    let dot = xs
        .iter()
        .zip(ys.iter())
        .fold(0.0, |acc, (x, y)| acc + x * y);
    Some(dot)
}

/// https://en.wikipedia.org/wiki/Sample_mean_and_covariance
pub fn covariance(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.len() != ys.len() || xs.len() <= 1 {
        return None;
    }

    let x_mean = mean(xs)?;
    let y_mean = mean(ys)?;

    let x_err: Vec<f64> = xs.iter().map(|x| x - x_mean).collect();
    let y_err: Vec<f64> = ys.iter().map(|y| y - y_mean).collect();

    let dot = dot(&x_err, &y_err)?;
    Some(dot / (xs.len() - 1) as f64)
}

 */

#[cfg(test)]
mod test {
    use crate::{Mean, Sum};

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

    // #[test]
    // fn variance() {
    //     assert_eq!(super::variance(&Vec::with_capacity(0)), None);

    //     let xs = vec![2.0, 2.0, 2.0, 2.0, 2.0];
    //     assert_eq!(super::variance(&xs), Some(0.0));

    //     let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    //     assert_eq!(super::variance(&xs), Some(2.0));

    //     let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    //     // assert_eq!(super::variance(&xs), Some(2.));
    //     // assert_eq!(super::variance(&xs), super::covariance(&xs, &xs)); // TODO: rescale
    // }

    // #[test]
    // fn dot() {
    //     let xs = vec![1.0];
    //     let ys = vec![];
    //     assert_eq!(super::dot(&xs, &ys), None);

    //     let xs = vec![1.0];
    //     let ys = vec![4.0, 5.0];
    //     assert_eq!(super::dot(&xs, &ys), None);

    //     let xs = vec![1.0, 2.0, 3.0];
    //     let ys = vec![4.0, 5.0, 6.0];
    //     assert_eq!(super::dot(&xs, &ys), Some(32.0));

    //     let xs = vec![1.0, 2.0, 3.0, 4.0];
    //     let ys = vec![4.0, 5.0, 6.0, 7.0];
    //     assert_eq!(super::dot(&xs, &ys), Some(60.0));

    //     let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    //     let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0];
    //     assert_eq!(super::dot(&xs, &ys), Some(100.0));

    //     let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    //     let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    //     assert_eq!(super::dot(&xs, &ys), Some(154.0));
    // }

    // #[test]
    // fn covariance() {
    //     let xs = vec![1.0];
    //     let ys = vec![];
    //     assert_eq!(super::covariance(&xs, &ys), None);

    //     let xs = vec![1.0];
    //     let ys = vec![4.0, 5.0];
    //     assert_eq!(super::covariance(&xs, &ys), None);

    //     let xs = vec![1.0, 2.0, 3.0];
    //     let ys = vec![4.0, 5.0, 6.0];
    //     assert_eq!(super::covariance(&xs, &ys), Some(1.0));

    //     // let xs = vec![1.0, 2.0, 3.0, 4.0];
    //     // let ys = vec![4.0, 5.0, 6.0, 7.0];
    //     // assert_eq!(super::covariance(&xs, &ys), Some(60.0));

    //     let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    //     let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0];
    //     assert_eq!(super::covariance(&xs, &ys), Some(2.5));

    //     // let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    //     // let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    //     // assert_eq!(super::covariance(&xs, &ys), Some(154.0));
    // }
}
