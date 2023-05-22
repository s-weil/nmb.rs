mod euler;
mod runge_kutta;

pub use euler::EulerSolver;
pub use runge_kutta::{Rk2Solver, Rk4Solver};

use nmbrs_algebra::VectorSpace;

#[derive(Debug, Clone, PartialEq)]
pub struct OdeState1D {
    pub t: f64,
    pub y: f64,
}
// TODO: add convenience methods and wrap it

pub trait OdeStepSolver1D {
    fn solve_step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64;
}

/// Solver for a first order ordninary differential equation ... // TODO
/// ```rust
/// use nmbrs_optimization::ode_solvers::{OdeSolver1D, OdeState1D};
///
/// // the initial value problem...
/// let f = |s: &OdeState1D| s.t.sin() * s.y ;
/// let initial_state = OdeState1D { t: 0.0, y: -1.0 };
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
/// let approximation_euler = EulerSolver.integrate(f, initial_state.clone(), t_end, n);
///
/// // ...has a convergence order of only 1
/// let upper_bound = 20.0 * dt.powi(1);
/// for s in approximation_euler {
///     assert!((sol(s.t) - s.y).abs() <= upper_bound);
/// }
///
/// // RungeKutta of order 2...
/// use nmbrs_optimization::ode_solvers::Rk2Solver;
/// let approximation_rk2 = Rk2Solver.integrate(f, initial_state.clone(), t_end, n);
///
/// // ...has a convergence order of 2
/// let upper_bound = 20.0 * dt.powi(2);
/// for s in approximation_rk2 {
///     assert!((sol(s.t) - s.y).abs() <= upper_bound);
/// }
///
/// // RungeKutta of order 4...
/// use nmbrs_optimization::ode_solvers::Rk4Solver;
/// let approximation_rk4 = Rk4Solver.integrate(f, initial_state, t_end, n);
///
/// // ...has a convergence order of 4
/// let upper_bound = 20.0 * dt.powi(4);
/// for s in approximation_rk4 {
///     assert!((sol(s.t) - s.y).abs() <= upper_bound);
/// }
/// ```
pub trait OdeSolver1D {
    fn integrate<F>(
        &self,
        f: F,
        initial_state: OdeState1D,
        t_end: f64,
        n: usize,
    ) -> Vec<OdeState1D>
    where
        F: Fn(&OdeState1D) -> f64;
}

impl<T> OdeSolver1D for T
where
    T: OdeStepSolver1D,
{
    fn integrate<F>(&self, f: F, initial_state: OdeState1D, t_end: f64, n: usize) -> Vec<OdeState1D>
    where
        F: Fn(&OdeState1D) -> f64,
    {
        integrate(self, f, initial_state, t_end, n)
    }
}

pub fn integrate<S, F>(
    solver: &S,
    f: F,
    initial_state: OdeState1D,
    t_end: f64,
    n: usize,
) -> Vec<OdeState1D>
where
    S: OdeStepSolver1D,
    F: Fn(&OdeState1D) -> f64,
{
    if t_end < initial_state.t || n < 1 {
        return Vec::with_capacity(0);
    }

    let dt = (t_end - initial_state.t) / n as f64;
    let mut ys = Vec::with_capacity(n + 1);
    ys.push(initial_state);

    for _ in 0..n {
        if let Some(state) = ys.last() {
            if state.t < t_end {
                let next_state = solver.solve_step(&f, state, dt);
                ys.push(next_state);
            }
        }
    }

    ys
}

pub trait Ode<V>: Fn(&TimeState<V>) -> V
where
    V: VectorSpace<f64>,
{
}

impl<F, V> Ode<V> for F
where
    F: Fn(&TimeState<V>) -> V,
    V: VectorSpace<f64>,
{
}

#[derive(Debug, Clone, PartialEq)]
pub struct TimeState<V>
where
    V: VectorSpace<f64>,
{
    pub t: f64,
    pub y: V,
}
// TODO: add convenience methods and wrap it

pub trait OdeStepSolver {
    fn solve_step<F, V>(
        &self,
        // for simplicity we assume that the domain and image of f is both V
        f: F,
        state: &TimeState<V>,
        dt: f64,
    ) -> TimeState<V>
    where
        F: Ode<V>,
        V: VectorSpace<f64>;
}

pub trait OdeSolver {
    // type VectorSpace: VectorSpace<f64>;

    fn integrate<F, V>(
        &self,
        f: F,
        initial_state: TimeState<V>,
        t_end: f64,
        n: usize,
    ) -> Vec<TimeState<V>>
    where
        F: Ode<V>,
        V: VectorSpace<f64>;
}

impl<T> OdeSolver for T
where
    T: OdeStepSolver,
{
    // type VectorSpace = V;

    fn integrate<F, V>(
        &self,
        f: F,
        initial_state: TimeState<V>,
        t_end: f64,
        n: usize,
    ) -> Vec<TimeState<V>>
    where
        F: Ode<V>,
        V: VectorSpace<f64>,
    {
        integrate2(self, f, initial_state, t_end, n)
    }
}

pub fn integrate2<S, F, V>(
    solver: &S,
    f: F,
    initial_state: TimeState<V>,
    t_end: f64,
    n: usize,
) -> Vec<TimeState<V>>
where
    S: OdeStepSolver,
    F: Ode<V>,
    V: VectorSpace<f64>,
{
    if t_end < initial_state.t || n < 1 {
        return Vec::with_capacity(0);
    }

    let dt = (t_end - initial_state.t) / n as f64;
    let mut ys = Vec::with_capacity(n + 1);
    ys.push(initial_state);

    for _ in 0..n {
        if let Some(state) = ys.last() {
            if state.t < t_end {
                let next_state = solver.solve_step(&f, state, dt);
                ys.push(next_state);
            }
        }
    }

    ys
}
