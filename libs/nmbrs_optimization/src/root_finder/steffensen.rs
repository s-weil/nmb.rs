use super::RootFinderConfig;

/*
PYTHON

def g(f: Func, x: float, fx: float) -> Func:
    """First-order divided difference function.

    Arguments:
        f: Function input to g
        x: Point at which to evaluate g
        fx: Function f evaluated at x
    """
    return lambda x: f(x + fx) / fx - 1

def steff(f: Func, x: float) -> Iterator[float]:
    """Steffenson algorithm for finding roots.

    This recursive generator yields the x_{n+1} value first then, when the generator iterates,
    it yields x_{n+2} from the next level of recursion.

    Arguments:
        f: Function whose root we are searching for
        x: Starting value upon first call, each level n that the function recurses x is x_n
    """
    while True:
        fx = f(x)
        gx = g(f, x, fx)(x)
        if gx == 0:
            break
        else:
            x = x - fx / gx    # Update to x_{n+1}
            yield x            # Yield value
 */

/// [Steffensen's method](https://en.wikipedia.org/wiki/Secant_method) for finding a root of a function `f`
/// is similiar to Newton's method, but uses a first-order divided difference function as approximation for the
/// derivative of `f`.
pub fn steffensen<F>(f: F, x0: f64, config: Option<RootFinderConfig>) -> Option<f64>
where
    F: Fn(f64) -> f64,
{
    let config = config.unwrap_or_default();
    let tol = config.tolerance;
    let max_iterations = config.max_iterations;

    let mut n_iterations = 0;
    let mut x = x0;

    while n_iterations < max_iterations {
        let f_x = f(x);

        if f_x.abs() < tol {
            return Some(x);
        }

        let df_x = f(x + f_x) / f_x - 1.0;

        if df_x.abs() < tol {
            return None;
        }

        let delta = -f_x / df_x;
        x += delta;

        if delta.abs() < tol {
            return Some(x);
        }

        n_iterations += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use std::f64::consts::SQRT_2;

    #[test]
    fn steffensen_root_quadratic() {
        let f = |x: f64| x * x - 2.0;

        // variant 1: start above the right root
        let root = super::steffensen(f, 3.0, None);
        assert_abs_diff_eq!(root.unwrap(), SQRT_2, epsilon = 1e-15);

        // variant 2: start below the right root
        let root = super::steffensen(f, 0.5, None);
        assert_abs_diff_eq!(root.unwrap(), -SQRT_2, epsilon = 1e-15);

        // variant 3: start above the left root
        let root = super::steffensen(f, -0.5, None);
        assert_abs_diff_eq!(root.unwrap(), -SQRT_2, epsilon = 1e-15);

        // variant 4: start in the middle between both roots
        let root = super::steffensen(f, 0.0, None);
        assert_abs_diff_eq!(root.unwrap(), -SQRT_2, epsilon = 1e-15);

        // variant 5: start left of the left root but sufficiently close
        let root = super::steffensen(f, -1.45, None);
        assert_abs_diff_eq!(root.unwrap(), -SQRT_2, epsilon = 1e-15);
    }

    #[test]
    fn steffenson_no_root() {
        let f = |x: f64| x * x - 2.0;

        // start left of the left root, but too far away
        let root = super::steffensen(f, -3.0, None);
        assert!(root.is_none());
    }
}
