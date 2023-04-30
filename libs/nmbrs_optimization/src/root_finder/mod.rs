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
        self.max_iterations = max_iterations;
        self
    }

    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
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

// pub trait RootFinder {
//     fn find_root<F>(&self, f: F, a: f64, b: f64, config: Option<RootFinderConfig>) -> f64
//     where
//         F: Fn(f64) -> f64;
// }
