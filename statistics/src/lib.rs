pub fn sum(xs: &[f64]) -> Option<f64> {
    if xs.is_empty() {
        return None;
    }

    let sum = xs.iter().fold(0.0, |acc, x| acc + x);
    Some(sum)
}

pub fn mean(xs: &[f64]) -> Option<f64> {
    let sum = sum(xs)?;
    Some(sum / xs.len() as f64)
}

/// NOTE: The variance is covered by the `Covariance` but provided as a more performant function.
pub fn variance(xs: &[f64]) -> Option<f64> {
    let mean = mean(xs)?;
    let mse = xs.iter().fold(0.0, |err, x| err + (x - mean).powi(2));
    Some(mse / xs.len() as f64)
}

#[cfg(test)]
mod test {

    #[test]
    fn sum() {
        assert_eq!(super::sum(&Vec::with_capacity(0)), None);
        assert_eq!(super::sum(&vec![]), None);

        let xs = vec![1.0, 1.0, 2.0];
        assert_eq!(super::sum(&xs), Some(4.0));

        let xs = vec![1.0, 2.0, 3.5];
        assert_eq!(super::sum(&xs), Some(6.5));
    }

    #[test]
    fn mean() {
        assert_eq!(super::mean(&Vec::with_capacity(0)), None);

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(super::mean(&xs), Some(3.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(super::mean(&xs), Some(3.5));
    }

    #[test]
    fn variance() {
        assert_eq!(super::variance(&Vec::with_capacity(0)), None);

        let xs = vec![2.0, 2.0, 2.0, 2.0, 2.0];
        assert_eq!(super::variance(&xs), Some(0.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(super::variance(&xs), Some(2.0));

        let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(super::variance(&xs), Some(3.5));
    }
}
