use super::RootFinderConfig;

/// implement the Newton method for root finding
/// https://en.wikipedia.org/wiki/Newton%27s_method
/// https://mathworld.wolfram.com/NewtonsMethod.html
pub fn newton<F, DF>(f: F, df: DF, x0: f64, config: Option<RootFinderConfig>) -> Option<f64>
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    let config = config.unwrap_or_default();
    let tol = config.tolerance.unwrap_or(1e-15);
    let max_iterations = config.max_iterations.unwrap_or(100);

    let mut x = x0;
    let mut df_x = df(x);

    // TODO: improve on thresholds, validations and error handling
    if df_x.abs() < 1e-15_f64.min(tol) {
        return None;
    }

    let mut f_x = f(x);

    let mut delta = -f_x / df_x;
    let mut n_iterations = 0;

    while delta.abs() > tol && n_iterations < max_iterations {
        x += delta;
        f_x = f(x);
        df_x = df(x);

        if df_x.abs() < 1e-15_f64.min(tol) {
            return None;
        }
        delta = -f_x / df_x;

        n_iterations += 1;
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_newton() {
        let f = |x: f64| x * x - 2.0;
        let df = |x: f64| 2.0 * x;
        let root = newton(f, df, 1.0, None);

        assert_abs_diff_eq!(root.unwrap(), 1.414213562373095, epsilon = 1e-15);
    }
}
