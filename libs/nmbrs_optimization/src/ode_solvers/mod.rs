mod runge_kutta;

pub trait InitialValueOde<const D: usize> {
    /// The function `f` in the problem `dy/dt = f(t, y(t))`
    fn f(t: f64, y: &[f64; D]) -> [f64; D];

    /// The intial value `y_0`.
    fn initial_value() -> [f64; D];
}
