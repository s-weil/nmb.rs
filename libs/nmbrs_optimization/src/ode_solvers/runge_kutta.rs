use crate::ode_solvers::{OdeState1D, OdeStepSolver1D};

// TODO: add D=1 case for scalar like odes;
// split in steps, maybe even have a trait for just solving one step

// see https://people.math.ethz.ch/~hiptmair/tmp/NUMODE11.pdf

// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

// todo: relax space to any vector space

pub struct Rk2Solver;

impl Rk2Solver {
    pub fn step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64,
    {
        // evaluate f at t, y_t
        let k1 = f(state);
        // approximate f(t + dt, y1) via Euler f(t + dt, y + dt * f(t , y))
        let k2 = f(&OdeState1D {
            t: state.t + dt,
            y: state.y + dt * k1,
        });
        let weighted_slope = (k1 + k2) / 2.0;

        // approximate y1 via Euler but the slope at t replaced by the mean of the slopes at t and t+dt, that is with k1 and k2

        OdeState1D {
            t: state.t + dt,
            y: state.y + dt * weighted_slope,
        }
    }
}

impl OdeStepSolver1D for Rk2Solver {
    fn solve_step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64,
    {
        self.step(f, state, dt)
    }
}

pub struct Rk4Solver;

impl Rk4Solver {
    pub fn step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64,
    {
        let k1 = f(state);
        let k2 = f(&OdeState1D {
            t: state.t + dt / 2.0,
            y: state.y + dt / 2.0 * k1,
        });
        let k3 = f(&OdeState1D {
            t: state.t + dt / 2.0,
            y: state.y + dt / 2.0 * k2,
        });
        let k4 = f(&OdeState1D {
            t: state.t + dt,
            y: state.y + dt * k3,
        });
        let weighted_slope = (k1 + 2.0 * (k2 + k3) + k4) / 6.0;

        // "Euler step"

        OdeState1D {
            t: state.t + dt,
            y: state.y + dt * weighted_slope,
        }
    }
}

impl OdeStepSolver1D for Rk4Solver {
    fn solve_step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64,
    {
        self.step(f, state, dt)
    }
}

// fn scalar<const D: usize>(scalar: f64, vec: &[f64; D]) -> [f64; D] {
//     let mut v = *vec;
//     for i in 0..D {
//         v[i] *= scalar;
//     }
//     v
// }

// fn add<const D: usize>(vec1: &[f64; D], vec2: &[f64; D]) -> [f64; D] {
//     let mut v = *vec1;
//     for i in 0..D {
//         v[i] += vec2[i];
//     }
//     v
// }

// pub fn runge_kutta_second_order<const D: usize, F>(
//     f: F,
//     y0: [f64; D],
//     t_end: f64,
//     n: usize,
// ) -> Vec<[f64; D]>
// where
//     F: Fn(T, &[T; D]) -> [T; D],
// {
//     let dt = t_end / n as f64;
//     let mut t = 0.0;

//     let mut ys = Vec::with_capacity(n + 1);

//     let mut y0 = y0;

//     while t < t_end {
//         let k1 = f(t, &y0);
//         let k2 = f(t + dt, &add(&y0, &scalar(dt, &k1)));
//         let y = add(&y0, &scalar(dt * 0.5, &(add(&k1, &k2))));
//         y0 = y;
//         t += dt;
//         ys.push(y);
//     }

//     assert_eq!(ys.len(), n);
//     ys
// }

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

    use crate::ode_solvers::OdeSolver1D;
    use crate::ode_solvers::OdeState1D;

    #[test]
    fn runge_kutta_second_order_1d_convegence() {
        // initial value problem
        let f = |s: &OdeState1D| s.y * s.t.sin();
        let y0 = -1.0;
        let initial_state = OdeState1D { t: 0.0, y: y0 };

        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let t_end = 10.0;

        for k in 5..15 {
            let n = 2_usize.pow(k);
            let ys = super::Rk2Solver.integrate(f, initial_state.clone(), t_end, n);

            let h = t_end / n as f64;
            let upper_bound = 5.0 * h.powi(2);

            for i in 0..n {
                let s_i = &ys[i];
                let sol_i = sol(s_i.t);
                let err_i = (sol_i - s_i.y).abs();
                assert!(err_i <= upper_bound);
            }
        }
    }

    #[test]
    fn runge_kutta_fourth_order_1d_convegence() {
        // initial value problem
        let f = |s: &OdeState1D| s.y * s.t.sin();
        let y0 = -1.0;
        let initial_state = OdeState1D { t: 0.0, y: y0 };
        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let t_end = 10.0;

        for k in 5..15 {
            let n = 2_usize.pow(k);
            let ys = super::Rk4Solver.integrate(f, initial_state.clone(), t_end, n);

            let h = t_end / n as f64;
            let upper_bound = 5.0 * h.powi(4);

            for i in 0..n {
                let s_i = &ys[i];
                let sol_i = sol(s_i.t);
                let err_i = (sol_i - s_i.y).abs();
                assert!(err_i <= upper_bound);
            }
        }
    }
}
