mod runge_kutta;

pub trait OdeInitialValueProblem<const N: usize> {
    /// The function `f` in the problem `dy/dt = f(t, y(t))`
    fn f(t: f64, y: &[f64; N]) -> [f64; N];

    fn initial_value();
}
