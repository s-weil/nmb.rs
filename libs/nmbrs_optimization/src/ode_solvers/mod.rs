mod euler;
mod runge_kutta;

use std::fmt::Display;

pub use euler::EulerSolver;
use nmbrs_algebra::VectorSpace;
pub use runge_kutta::{Rk2Solver, Rk4Solver};

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

// pub trait OdeSystem {
//     type Space: VectorSpace;

//     fn ode(time_state: &TimeState<Self::Space>) -> Self::Space;
// }

// impl<F, V: VectorSpace> OdeSystem for F
// where
//     F: Fn(&TimeState<V>) -> V,
// {
//     type Space = V;
//     fn ode(time_state: &TimeState<Self::Space>) -> Self::Space {
//         Self(time_state)
//     }
// }

// impl<F, V: VectorSpace> Ode for F
// where
//     F: Fn(&TimeState<V>) -> V,
// {
//     type Space = V;
// }

// for simplicity we assume that the domain and image of f is both V
// pub trait OdeSystem: Fn(&TimeState<Self::Space>) -> Self::Space {
//     type Space: VectorSpace;
// }

// impl<F, V> OdeSystem for F
// where
//     F: Fn(&TimeState<V>) -> V,
//     V: VectorSpace,
// {
//     type Space = V;
// }

// for simplicity we assume that the domain and image of f is both V
pub trait OdeSystem2<V>: Fn(&TimeState<V>) -> V
where
    V: VectorSpace,
{
}

impl<F, V> OdeSystem2<V> for F
where
    F: Fn(&TimeState<V>) -> V,
    V: VectorSpace,
{
}

// TODO: which version is nicer/more convenient?

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
    V: std::fmt::Debug,
    V::Field: std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{:?}: {:?}]", self.t, self.y)
    }
    // fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    //     write!(f, "({}: {})", self.t, self.y)
    // }
}

// pub trait OdeStepSolver {
//     fn solve_step<S>(
//         &self,
//         // for simplicity we assume that the domain and image of f is both V
//         f: &S,
//         state: &TimeState<S::Space>,
//         dt: <S::Space as VectorSpace>::Field,
//     ) -> TimeState<S::Space>
//     where
//         S: OdeSystem,
//         S::Space: Clone,
//         <S::Space as VectorSpace>::Field: Clone;
// }

pub trait OdeStepSolver2<S, V>
where
    S: OdeSystem2<V>,
    V: VectorSpace + Clone,
    V::Field: Clone,
{
    fn solve_step(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V>;
}

// pub trait OdeSolver {
//     // type VectorSpace: VectorSpace<f64>;

//     fn integrate<S>(
//         &self,
//         f: &S,
//         initial_state: TimeState<S::Space>,
//         t_end: <S::Space as VectorSpace>::Field,
//         n: usize,
//     ) -> Vec<TimeState<S::Space>>
//     where
//         S: OdeSystem,
//         S::Space: Clone,
//         <S::Space as VectorSpace>::Field: Clone + From<i8> + PartialOrd;
// }

pub trait OdeSolver2<S, V>
where
    S: OdeSystem2<V>,
    V: VectorSpace + Clone,
    V::Field: Clone,
{
    fn integrate(
        &self,
        f: &S,
        initial_state: TimeState<V>,
        t_end: V::Field,
        n: usize,
    ) -> Vec<TimeState<V>>;
}

// impl<T> OdeSolver for T
// where
//     T: OdeStepSolver,
// {
//     // type VectorSpace = V;

//     fn integrate<S>(
//         &self,
//         f: &S,
//         initial_state: TimeState<S::Space>,
//         t_end: <S::Space as VectorSpace>::Field,
//         n: usize,
//     ) -> Vec<TimeState<S::Space>>
//     where
//         S: OdeSystem,
//         S::Space: Clone,
//         <S::Space as VectorSpace>::Field: Clone + From<i8> + PartialOrd,
//     {
//         integrate2(self, f, initial_state, t_end, n)
//     }
// }

impl<T, S, V> OdeSolver2<S, V> for T
where
    T: OdeStepSolver2<S, V>,
    S: OdeSystem2<V>,
    V: VectorSpace + Clone,
    V::Field: Clone + From<i8> + PartialOrd,
{
    fn integrate(
        &self,
        f: &S,
        initial_state: TimeState<V>,
        t_end: V::Field,
        n: usize,
    ) -> Vec<TimeState<V>> {
        integrate3(self, f, initial_state, t_end, n)
    }
}

// pub fn integrate2<X, S>(
//     solver: &X,
//     f: &S,
//     initial_state: TimeState<S::Space>,
//     t_end: <S::Space as VectorSpace>::Field,
//     n: usize,
// ) -> Vec<TimeState<S::Space>>
// where
//     X: OdeStepSolver,
//     S: OdeSystem,
//     S::Space: Clone,
//     <S::Space as VectorSpace>::Field: Clone + From<i8> + PartialOrd,
// {
//     if t_end < initial_state.t || n < 1 {
//         return Vec::with_capacity(0);
//     }

//     let dt = (t_end.clone() - initial_state.t.clone()) / (n as i8).into();
//     let mut ys = Vec::with_capacity(n + 1);
//     ys.push(initial_state);

//     for _ in 0..n {
//         if let Some(state) = ys.last() {
//             if state.t < t_end {
//                 let next_state = solver.solve_step(f, state, dt.clone());
//                 ys.push(next_state);
//             }
//         }
//     }

//     ys
// }

pub fn integrate3<X, S, V>(
    solver: &X,
    f: &S,
    initial_state: TimeState<V>,
    t_end: V::Field,
    n: usize,
) -> Vec<TimeState<V>>
where
    X: OdeStepSolver2<S, V>,
    S: OdeSystem2<V>,
    V: VectorSpace + Clone,
    V::Field: Clone + From<i8> + PartialOrd,
{
    if t_end < initial_state.t || n < 1 {
        return Vec::with_capacity(0);
    }

    let dt = (t_end.clone() - initial_state.t.clone()) / (n as i8).into();
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
