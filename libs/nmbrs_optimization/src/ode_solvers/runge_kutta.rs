use super::{OdeStepSolver, OdeSystem, TimeState};
use nmbrs_algebra::{NumericField, VectorSpace};

// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

fn weight<F>(denominator: usize) -> F
where
    F: NumericField + From<i32>,
{
    F::one() / F::from(denominator as i32)
}

/// The [Runge Kutta Method](https://en.wikipedia.org/wiki/Runge-Kutta_methods)
/// of order 2.
pub struct Rk2Solver;

impl Rk2Solver {
    pub fn step<S, V>(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V>
    where
        S: OdeSystem<V>,
        V: VectorSpace + Clone,
        V::Field: Clone + From<i32>,
    {
        // in short:
        // let k1 = f(state);
        // let k2 = f(&OdeState1D {
        //     t: state.t + dt,
        //     y: state.y + dt * k1,
        // });
        // let weighted_slope = (k1 + k2) / 2.0;

        // k1: evaluate f at t, y_t
        let k1 = f(state);

        // k2: approximate f(t + dt, y1) via Euler f(t + dt, y + dt * f(t , y))
        let t_step = state.t.clone() + dt.clone();
        let k2 = f(&TimeState {
            t: t_step.clone(),
            y: state.y.clone() + f(state) * dt.clone(),
        });

        // approximate y1 via Euler but the slope at t replaced by the mean of the slopes at t and t+dt,
        // that is with the average of k1 and k2
        let weighted_slope = (k1 + k2) * weight(2);

        TimeState {
            t: t_step,
            y: state.y.clone() + weighted_slope * dt,
        }
    }
}

impl<S, V> OdeStepSolver<S, V> for Rk2Solver
where
    S: OdeSystem<V>,
    V: VectorSpace + Clone,
    V::Field: Clone + From<i32>,
{
    fn solve_step(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V> {
        self.step(f, state, dt)
    }
}

/// The [Runge Kutta Method](https://en.wikipedia.org/wiki/Runge-Kutta_methods)
/// of order 4.
pub struct Rk4Solver;

impl Rk4Solver {
    pub fn step<S, V>(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V>
    where
        S: OdeSystem<V>,
        V: VectorSpace + Clone,
        V::Field: Clone + From<i32>,
    {
        // in short:
        // let k1 = f(state);
        // let k2 = f(&TimeState {
        //     t: state.t + dt / 2.0,
        //     y: state.y + dt / 2.0 * k1,
        // });
        // let k3 = f(&TimeState {
        //     t: state.t + dt / 2.0,
        //     y: state.y + dt / 2.0 * k2,
        // });
        // let k4 = f(&TimeState {
        //     t: state.t + dt,
        //     y: state.y + dt * k3,
        // });
        // let weighted_slope = (k1 + 2.0 * (k2 + k3) + k4) / 6.0;

        let k1 = f(state);

        // k2 & k3: take approximate derivatives at t + dt/2
        let dt_mid = dt.clone() * weight(2);
        let t_mid = state.t.clone() + dt_mid.clone();
        let k2: V = f(&TimeState {
            t: t_mid.clone(),
            y: state.y.clone() + k1.clone() * dt_mid.clone(),
        });
        let k3 = f(&TimeState {
            t: t_mid,
            y: state.y.clone() + k2.clone() * dt_mid,
        });
        let k_mid = (k2 + k3.clone()) * V::Field::from(2);

        let t_step = state.t.clone() + dt.clone();
        let k4 = f(&TimeState {
            t: t_step.clone(),
            y: state.y.clone() + k3 * dt.clone(),
        });

        let weighted_slope = (k1 + k_mid + k4) * weight(6);

        // "Euler step"
        TimeState {
            t: t_step,
            y: state.y.clone() + weighted_slope * dt,
        }
    }
}

impl<S, V> OdeStepSolver<S, V> for Rk4Solver
where
    S: OdeSystem<V>,
    V: VectorSpace + Clone,
    V::Field: Clone + From<i32>,
{
    fn solve_step(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V> {
        self.step(f, state, dt)
    }
}

#[cfg(test)]
mod tests {

    // #[test]
    // fn runge_kutta_second_order_1d_convegence_order() {
    //     // initial value problem
    //     let a = 1.0;
    //     let b = 0.0;
    //     let f = |t, y| a * t * y + b;
    //     let y0 = 0.5;

    //     // solution
    //     let sol = |t: f64| 0.5 * (t * t).exp();

    //     let t_end = 2.0;

    //     for k in 5..10 {
    //         let n = 2_usize.pow(k);
    //         let ys = super::runge_kutta_second_order_1d(f, y0, t_end, n);

    //         let h = t_end / n as f64;
    //         let upper_bound = 100.0 * h.powi(2);

    //         for i in 0..n {
    //             let t = i as f64 * h;
    //             let sol_i = sol(t);
    //             let (ti, yi) = ys[i];
    //             let err_i = (sol_i - yi).abs();
    //             dbg!(err_i, t, ti, sol_i, ys[i], i, n);
    //             assert!(err_i <= upper_bound);
    //         }
    //     }
    // }

    use crate::ode_solvers::{OdeSolver, TimeState};

    #[test]
    fn runge_kutta_second_order_1d_convegence() {
        // initial value problem
        let f = |s: &TimeState<f64>| s.y * s.t.sin();
        let y0 = -1.0;
        let initial_state: TimeState<f64> = TimeState { t: 0.0, y: y0 };

        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let t_end = 10.0;

        for k in 5..15 {
            let n = 2_usize.pow(k);
            let ys = super::Rk2Solver.integrate(&f, initial_state.clone(), t_end, n);

            let dt: f64 = t_end / n as f64;
            let upper_bound = 5.0 * dt.powi(2);

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

    #[test]
    fn runge_kutta_fourth_order_1d_convegence() {
        // initial value problem
        let f = |s: &TimeState<f64>| s.y * s.t.sin();
        let y0 = -1.0;
        let initial_state = TimeState { t: 0.0, y: y0 };

        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let t_end = 10.0;

        for k in 5..15 {
            let n = 2_usize.pow(k);
            let ys = super::Rk4Solver.integrate(&f, initial_state.clone(), t_end, n);

            let dt: f64 = t_end / n as f64;
            let upper_bound = 5.0 * dt.powi(4);

            for i in 0..n {
                let s_i = &ys[i];
                let sol_i = sol(s_i.t);
                let err_i = (sol_i - s_i.y).abs();
                assert!(
                    err_i <= upper_bound,
                    "error {} exceeded threshold {} ({})",
                    err_i,
                    upper_bound,
                    n
                );
            }
        }
    }
}
