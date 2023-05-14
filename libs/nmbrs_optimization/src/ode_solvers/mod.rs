mod runge_kutta;

type T = f64;

// pub trait Ode<const D: usize> {
//     /// The function `f` in the problem `dy/dt = f(t, y)`
//     fn f(t: T, y: &[T; D]) -> [T; D];
// }

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
    // TODO: validation of inputs
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
