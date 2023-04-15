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
    if xs.len() != ys.len() || xs.len() == 1 {
        return None;
    }

    let x_mean = mean(xs)?;
    let y_mean = mean(ys)?;

    let x_err: Vec<f64> = xs.iter().map(|x| x - x_mean).collect();
    let y_err: Vec<f64> = ys.iter().map(|y| y - y_mean).collect();

    let dot = dot(&x_err, &y_err)?;
    Some(dot / (xs.len() - 1) as f64)
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

        // let xs = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        // let ys = vec![4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        // assert_eq!(super::covariance(&xs, &ys), Some(154.0));
    }
}
