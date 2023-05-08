mod runge_kutta;

type T = f64;

pub trait Ode<const D: usize> {
    /// The function `f` in the problem `dy/dt = f(t, y)`
    fn f(t: T, y: &[T; D]) -> [T; D];
}

pub trait InitialValueOde<const D: usize> {
    /// The function `f` in the problem `dy/dt = f(t, y(t))`
    fn f(t: f64, y: &[f64; D]) -> [f64; D];

    /// The intial value `y_0`.
    fn y0() -> [f64; D];
}

pub trait OdeSolver<const D: usize> {
    fn solve<F>(ode: F, t: f64, n: usize) -> Vec<[f64; D]>
    where
        F: InitialValueOde<D>;
}
