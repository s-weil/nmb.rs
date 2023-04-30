mod bisection;
mod newton;

pub use bisection::bisection;
pub use newton::newton;

#[derive(Debug, Clone)]
pub struct RootFinderConfig {
    pub max_iterations: usize,
    pub tolerance: f64,
}

impl RootFinderConfig {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        if max_iterations == 0 {
            panic!("max_iterations must be greater than 0");
        }
        self.max_iterations = max_iterations;
        self
    }

    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        if tolerance <= 0.0 {
            panic!("tolerance must be greater than 0");
        }
        self.tolerance = tolerance;
        self
    }
}

impl Default for RootFinderConfig {
    fn default() -> Self {
        Self {
            max_iterations: (100),
            tolerance: (1e-15),
        }
    }
}

/// A root solver for finding a solution to the equation $f(x) = 0$.
/// See https://en.wikipedia.org/wiki/Root-finding_algorithms
///
/// ```rust
/// use nmbrs_optimization::root_finder::RootSolver;
/// // define the function $f(x) = x^2 - 2$ for which we want to find the root
/// let f = |x: f64| x * x - 2.0;
///
/// // use the Newton Raphson algorithm which requires the derivative of f and a guess for the starting point
/// use nmbrs_optimization::root_finder::DerivativeSolver;
/// let df = |x: f64| 2.0 * x;
/// let root = DerivativeSolver::newton_raphson(f, df, 3.0).try_find_root(None);
/// assert!( (root.unwrap() - 2.0_f64.sqrt()).abs() < 1e-15);
/// // if you start with a guess that is too far away from the root or at a point where $df=0$, the algorithm might fail
/// assert!(DerivativeSolver::newton_raphson(f, df, 0.0).try_find_root(None).is_none());
///
/// // use the bisection algorithm which requires a bracketing interval
/// use nmbrs_optimization::root_finder::BracketingSolver;
/// let root = BracketingSolver::bisection(f, 0.0, 3.0).try_find_root(None);
/// assert!( (root.unwrap() - 2.0_f64.sqrt()).abs() < 1e-15);
/// // if you select an interval for which both $f(a)$ and $f(b)$ have the same sign, the algorithm will fail
/// assert!(BracketingSolver::bisection(f, -1.0, 1.0).try_find_root(None).is_none());
/// ```
pub trait RootSolver {
    fn try_find_root(&self, config: Option<RootFinderConfig>) -> Option<f64>;
}

pub enum BracketingSolver<F> {
    Bisection { f: F, a: f64, b: f64 },
}

impl<F> BracketingSolver<F>
where
    F: Fn(f64) -> f64,
{
    pub fn bisection(f: F, a: f64, b: f64) -> Self {
        Self::Bisection { f, a, b }
    }
}

impl<F> RootSolver for BracketingSolver<F>
where
    F: Fn(f64) -> f64,
{
    fn try_find_root(&self, config: Option<RootFinderConfig>) -> Option<f64> {
        match self {
            Self::Bisection { f, a, b } => bisection(f, *a, *b, config),
        }
    }
}

pub enum DerivativeSolver<F, DF> {
    NewtonRaphson { f: F, df: DF, x0: f64 },
}

impl<F, DF> RootSolver for DerivativeSolver<F, DF>
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    fn try_find_root(&self, config: Option<RootFinderConfig>) -> Option<f64> {
        match self {
            Self::NewtonRaphson { f, df, x0 } => newton(f, df, *x0, config),
        }
    }
}

impl<F, DF> DerivativeSolver<F, DF>
where
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
{
    pub fn newton_raphson(f: F, df: DF, x0: f64) -> Self {
        Self::NewtonRaphson { f, df, x0 }
    }
}
