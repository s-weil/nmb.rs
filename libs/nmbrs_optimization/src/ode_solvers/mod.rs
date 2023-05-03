mod runge_kutta;

<<<<<<< HEAD
pub trait InitialValueOde<const D: usize> {
    /// The function `f` in the problem `dy/dt = f(t, y(t))`
    fn f(t: f64, y: &[f64; D]) -> [f64; D];

    /// The intial value `y_0`.
    fn initial_value() -> [f64; D];
=======
pub trait OdeInitialValueProblem<const N: usize> {
    /// The function `f` in the problem `dy/dt = f(t, y(t))`
    fn f(t: f64, y: &[f64; N]) -> [f64; N];

    fn initial_value();
>>>>>>> bf1a759a813b1c94dbf7a7e1065b2a08dc00cdfa
}
