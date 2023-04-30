mod bisection;
mod newton;

#[derive(Debug)]
pub struct RootFinderConfig {
    pub max_iterations: Option<usize>,
    pub tolerance: Option<f64>,
}

impl Default for RootFinderConfig {
    fn default() -> Self {
        Self {
            max_iterations: Some(100),
            tolerance: Some(1e-15),
        }
    }
}
