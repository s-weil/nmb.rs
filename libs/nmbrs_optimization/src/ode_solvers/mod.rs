mod euler;
mod runge_kutta;
pub use euler::EulerSolver;
use nmbrs_algebra::VectorSpace;
pub use runge_kutta::{Rk2Solver, Rk4Solver};
use std::fmt::{Debug, Display};

// for simplicity we assume that the domain and image of f is both V
pub trait OdeSystem<V>: Fn(&TimeState<V>) -> V
where
    V: VectorSpace,
{
}

impl<F, V> OdeSystem<V> for F
where
    F: Fn(&TimeState<V>) -> V,
    V: VectorSpace,
{
}

// #[derive(Debug, Clone, PartialEq)]
pub struct TimeState<V>
where
    V: VectorSpace,
{
    pub t: <V as VectorSpace>::Field,
    pub y: V,
}
// TODO: add convenience methods and wrap it

impl<V: VectorSpace> Clone for TimeState<V>
where
    V: Clone,
    V::Field: Clone,
{
    fn clone(&self) -> Self {
        TimeState::<V> {
            t: self.t.clone(),
            y: self.y.clone(),
        }
    }
}

impl<V: VectorSpace> Display for TimeState<V>
where
    V: Display,
    V::Field: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({}: {})", self.t, self.y)
    }
}

impl<V: VectorSpace> std::fmt::Debug for TimeState<V>
where
    V: Debug,
    V::Field: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{:?}: {:?}]", self.t, self.y)
    }
}

pub trait OdeStepSolver<S, V>
where
    S: OdeSystem<V>,
    V: VectorSpace + Clone,
    V::Field: Clone,
{
    fn solve_step(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V>;
}

/// [Numerical solver](https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations)
/// for a first order ordinary differential equation (ODE) with initial value (IVP).
/// ```rust
/// use nmbrs_optimization::ode_solvers::{OdeSolver, TimeState};
///
/// // the initial value problem...
/// let f = |s: &TimeState<f64>| s.y * s.t.sin();
/// let y0 = -1.0;
/// let initial_state: TimeState<f64> = TimeState { t: 0.0, y: y0 };
///
/// // ...has the solution
/// let sol = |t: f64| -(1.0 - t.cos()).exp();
///
/// let n = 1_000;
/// let t_end = 10.0;
/// let dt = t_end / n as f64;
///
///
/// // Euler of order 2...
/// use nmbrs_optimization::ode_solvers::EulerSolver;
/// let approximation_euler = EulerSolver.integrate(&f, initial_state.clone(), t_end, n);
///
/// // ...has a convergence order of only 1
/// let upper_bound = 20.0 * dt.powi(1);
/// for s in approximation_euler {
///     assert!((sol(s.t) - s.y).abs() <= upper_bound);
/// }
///
/// // RungeKutta of order 2...
/// use nmbrs_optimization::ode_solvers::Rk2Solver;
/// let approximation_rk2 = Rk2Solver.integrate(&f, initial_state.clone(), t_end, n);
///
/// // ...has a convergence order of 2
/// let upper_bound = 20.0 * dt.powi(2);
/// for s in approximation_rk2 {
///     assert!((sol(s.t) - s.y).abs() <= upper_bound);
/// }
///
/// // RungeKutta of order 4...
/// use nmbrs_optimization::ode_solvers::Rk4Solver;
/// let approximation_rk4 = Rk4Solver.integrate(&f, initial_state, t_end, n);
///
/// // ...has a convergence order of 4
/// let upper_bound = 20.0 * dt.powi(4);
/// for s in approximation_rk4 {
///     assert!((sol(s.t) - s.y).abs() <= upper_bound);
/// }
/// ```
pub trait OdeSolver<S, V>
where
    S: OdeSystem<V>,
    V: VectorSpace,
{
    fn integrate(
        &self,
        f: &S,
        initial_state: TimeState<V>,
        t_end: V::Field,
        n: usize,
    ) -> Vec<TimeState<V>>;
}

impl<T, S, V> OdeSolver<S, V> for T
where
    T: OdeStepSolver<S, V>,
    S: OdeSystem<V>,
    V: VectorSpace + Clone,
    V::Field: Clone + PartialOrd + From<i32>,
{
    fn integrate(
        &self,
        f: &S,
        initial_state: TimeState<V>,
        t_end: V::Field,
        n: usize,
    ) -> Vec<TimeState<V>> {
        integrate(self, f, initial_state, t_end, n)
    }
}

pub fn integrate<X, S, V>(
    solver: &X,
    f: &S,
    initial_state: TimeState<V>,
    t_end: V::Field,
    n: usize,
) -> Vec<TimeState<V>>
where
    X: OdeStepSolver<S, V>,
    S: OdeSystem<V>,
    V: VectorSpace + Clone,
    V::Field: Clone + PartialOrd + From<i32>,
{
    if t_end < initial_state.t || n < 1 {
        return Vec::with_capacity(0);
    }

    let dt = (t_end.clone() - initial_state.t.clone()) / (n as i32).into();
    let mut ys = Vec::with_capacity(n + 1);
    ys.push(initial_state);

    for _ in 0..n {
        if let Some(state) = ys.last() {
            if state.t < t_end {
                let next_state = solver.solve_step(f, state, dt.clone());
                ys.push(next_state);
            }
        }
    }

    ys
}
