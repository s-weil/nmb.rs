use super::RootFinderConfig;

/*
PYTHON

def secant_method(f, x0, x1, iterations):
    """Return the root calculated using the secant method."""
    for i in range(iterations):
        x2 = x1 - f(x1) * (x1 - x0) / float(f(x1) - f(x0))
        x0, x1 = x1, x2
        # Apply a stopping criterion here
    return x2
*/

/// Steffensen's method for finding a root of a function.
/// https://en.wikipedia.org/wiki/Steffensen%27s_method
pub fn secant<F>(f: F, x0: f64, x1: f64, config: Option<RootFinderConfig>) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let config = config.unwrap_or_default();
    let tol = config.tolerance;
    let max_iterations = config.max_iterations;

    let mut n_iterations = 0;
    let mut x0 = x0;
    let mut x1 = x1;

    if (x0 - x1).abs() < tol {
        panic!("initially guessed x0 and x1 are too close to each other");
        return None;
    }

    while n_iterations < max_iterations {
        let f_1 = f(x1);
        let x_diff = x1 - x0;

        if f_1.abs() < tol || x_diff.abs() < tol {
            return Some(x1);
        }

        let f_diff = f_1 - f(x0);

        if f_diff.abs() < tol {
            return None;
        }

        let x2 = x1 - f_1 * x_diff / f_diff;

        x0 = x1;
        x1 = x2;
        n_iterations += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use std::f64::consts::SQRT_2;

    #[test]
    fn secant_root_quadratic() {
        let f = |x: f64| x * x - 2.0;

        // variant 1: start above the right root
        let root = super::secant(f, 2.0, 4.0, None);
        assert_abs_diff_eq!(root.unwrap(), SQRT_2, epsilon = 1e-15);

        // variant 2: start below and above the right root
        let root = super::secant(f, 0.5, 3.0, None);
        assert_abs_diff_eq!(root.unwrap(), SQRT_2, epsilon = 1e-15);

        // variant 3: start in the middle between both roots
        let root = super::secant(f, 0.0, -1.0, None);
        assert_abs_diff_eq!(root.unwrap(), -SQRT_2, epsilon = 1e-15);

        // variant 4: start left of the left root
        let root = super::secant(f, -2., -0.4, None);
        assert_abs_diff_eq!(root.unwrap(), -SQRT_2, epsilon = 1e-15);
    }

    #[test]
    fn secant_no_root() {
        let f = |x: f64| x * x - 2.0;

        // guess symmetrically around the point with zero derivative
        let root = super::secant(f, -0.5, 0.5, None);
        assert!(root.is_none());

        // sguess symmetrically around the point with zero derivative
        let root = super::secant(f, -3.0, 3.0, None);
        assert!(root.is_none());
    }
}
