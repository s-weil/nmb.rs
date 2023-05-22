use std::marker::PhantomData;

use nmbrs_algebra::VectorSpace;

use crate::ode_solvers::{Ode, OdeStepSolver, TimeState};

pub struct EulerSolver<V>
where
    V: VectorSpace<f64>,
{
    v: PhantomData<V>,
}

// impl EulerSolver {
//     pub fn step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
//     where
//         F: Fn(&OdeState1D) -> f64,
//     {
//         OdeState1D {
//             t: state.t + dt,
//             y: state.y + dt * f(state),
//         }
//     }
// }

impl<V> EulerSolver<V>
where
    V: VectorSpace<f64> + Clone,
{
    pub fn step<F>(&self, f: F, state: &TimeState<V>, dt: f64) -> TimeState<V>
    where
        F: Ode<V>,
    {
        TimeState {
            t: state.t + dt,
            y: state.y.clone() + f(state) * dt,
        }
    }
}

impl<V> OdeStepSolver for EulerSolver<V>
where
    V: VectorSpace<f64> + Clone,
{
    fn solve_step<F, V>(&self, f: F, state: &TimeState<V>, dt: f64) -> TimeState<V>
    where
        F: Ode<V>,
    {
        self.step(f, state, dt)
    }
}

#[cfg(test)]
mod tests {

    use crate::ode_solvers::{OdeSolver1D, OdeState1D};

    #[test]
    fn euler_1d_convegence() {
        // initial value problem
        let f = |s: &OdeState1D| s.y * s.t.sin();
        let y0 = -1.0;
        let initial_state = OdeState1D { t: 0.0, y: y0 };
        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let t_end = 10.0;

        for k in 5..15 {
            let n = 2_usize.pow(k);
            let ys = super::EulerSolver.integrate(f, initial_state.clone(), t_end, n);

            let h = t_end / n as f64;
            let upper_bound = 20.0 * h;

            for i in 0..n {
                let s_i = &ys[i];
                let sol_i = sol(s_i.t);
                let err_i = (sol_i - s_i.y).abs();
                assert!(err_i <= upper_bound);
            }
        }
    }
}
