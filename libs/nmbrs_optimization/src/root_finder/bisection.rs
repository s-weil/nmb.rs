

use crate::root_finder::RootFinderConfig;

/// implementation of the bisection algorithm
/// https://en.wikipedia.org/wiki/Bisection_method
/// https://mathworld.wolfram.com/Bisection.htm
/// https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics/RootFinding/Bisection.csl
pub fn bisection<F>(f: F, a: f64, b: f64, config: Option<RootFinderConfig>) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let mut a = a;
    let mut b = b;

    if b < a {
        // std::mem::swap(&mut a, &mut b);
        return None;
    }

    let mut f_a = f(a);
    let mut f_b = f(b);

    if f_a * f_b > 0.0 {
        // TODO: proper error handling
        // panic!("f(a) and f(b) must have opposite signs");
        return None;
    }

    let config = config.unwrap_or_default();

    let tol = config.tolerance;
    let max_iterations = config.max_iterations;
    // .max_iterations
    // .unwrap_or(tol.log2().ceil() );

    if f_a.abs() < tol {
        return Some(a);
    }
    if f_b.abs() < tol {
        return Some(b);
    }

    let mut mid: f64 = (a + b) / 2.0;
    let mut f_mid = f(mid);
    let mut iterations = 0;

    let mut delta = b - a;
    while delta > tol && f_mid.abs() > tol && iterations < max_iterations {
        if f_a * f_mid < 0.0 {
            b = mid;
            f_b = f_mid;
        } else {
            a = mid;
            f_a = f_mid;
        }
        delta = b - a;
        mid = (a + b) / 2.0;
        f_mid = f(mid);
        iterations += 1;
    }
    Some(mid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_bisection_default() {
        let f = |x: f64| x * x - 2.0;

        // search for sqrt(2) in the interval [1, 2]
        let root = bisection(f, 1.0, 2.0, None);
        assert_abs_diff_eq!(root.unwrap(), 1.414213562373095, epsilon = 1e-15);

        // search for sqrt(2) in the interval [-2, 0]
        let root = bisection(f, -2.0, 0.0, None);
        assert_abs_diff_eq!(root.unwrap(), -1.414213562373095, epsilon = 1e-15);
    }

    #[test]
    fn test_bisection_no_root() {
        let f = |x: f64| x * x - 2.0;

        // no root in the interval: f(3) > 0 and f(4) > 0
        assert!(f(3.0) > 0.0);
        assert!(f(4.0) > 0.0);
        let root = bisection(f, 3.0, 4.0, None);
        assert!(root.is_none());

        // no root in the interval: f(-1) < 0 and f(1) < 0
        assert!(f(1.0) < 0.0);
        assert!(f(-1.0) < 0.0);
        let root = bisection(f, -1.0, -1.0, None);
        assert!(root.is_none());

        // cannot find root in the interval even though it exists: f(-2) > 0 and f(2) > 0
        assert!(f(2.0) > 0.0);
        assert!(f(2.0) > 0.0);
        let root = bisection(f, 2.0, 2.0, None);
        assert!(root.is_none());
        // TODO: provide version with randomized evaluations in the interval in order to find the root
    }
}
