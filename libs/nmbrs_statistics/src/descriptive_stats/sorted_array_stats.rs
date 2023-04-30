use crate::AsSlice;
use nmbrs_algebra::{MidPoint, NumericField};

/*
Sorted array statistics provides routines optimized for an array sorting ascendingly.
Especially order-statistics are very efficient this way, some even with constant time complexity.
Source: https://numerics.mathdotnet.com/DescriptiveStatistics
*/

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
mod tests {
    use super::Percentile;

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
}
