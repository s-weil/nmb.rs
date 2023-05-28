use crate::ode_solvers::{OdeStepSolver2, OdeSystem2, TimeState};
use nmbrs_algebra::VectorSpace;

pub struct EulerSolver;

impl EulerSolver {
    pub fn step<S, V>(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V>
    where
        S: OdeSystem2<V>,
        V: VectorSpace + Clone,
        V::Field: Clone,
    {
        TimeState {
            t: state.t.clone() + dt.clone(),
            y: state.y.clone() + f(state) * dt.clone(),
        }
    }
}

impl<S, V> OdeStepSolver2<S, V> for EulerSolver
where
    S: OdeSystem2<V>,
    V: VectorSpace + Clone + std::fmt::Debug,
    V::Field: Clone + std::fmt::Debug,
{
    fn solve_step(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V> {
        self.step(f, state, dt)
    }
}

#[cfg(test)]
mod tests {
    use crate::ode_solvers::{OdeSolver2, TimeState};

    #[test]
    fn euler_1d_convegence() {
        // initial value problem
        let f = |s: &TimeState<f64>| s.y * s.t.sin();
        let y0 = -1.0;
        let initial_state = TimeState { t: 0.0, y: y0 };
        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let solver = super::EulerSolver;
        let t_end = 10.0;

        for k in 5..15 {
            let n = 2_usize.pow(k);
            let ys: Vec<TimeState<f64>> = solver.integrate(&f, initial_state.clone(), t_end, n);

            let h = t_end / n as f64;
            let upper_bound = 20.0 * h;

            for i in 0..n {
                let s_i = &ys[i];
                let sol_i = sol(s_i.t);
                let err_i = (sol_i - s_i.y).abs();
                assert!(
                    err_i <= upper_bound,
                    "error {} exceeded threshold {}",
                    err_i,
                    upper_bound
                );
            }
        }
    }
}
