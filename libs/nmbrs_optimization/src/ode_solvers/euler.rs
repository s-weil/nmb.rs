use crate::ode_solvers::{OdeStepSolver2, OdeSystem2, TimeState};
use nmbrs_algebra::VectorSpace;
use std::marker::PhantomData;

// pub struct EulerSolver<S, V>
// where
//     S: OdeSystem<Space = V>,
//     V: VectorSpace,
// {
//     v: PhantomData<V>,
//     s: PhantomData<S>,
// }

// // impl EulerSolver {
// //     pub fn step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
// //     where
// //         F: Fn(&OdeState1D) -> f64,
// //     {
// //         OdeState1D {
// //             t: state.t + dt,
// //             y: state.y + dt * f(state),
// //         }
// //     }
// // }

// impl<S, V> EulerSolver<S, V>
// where
//     // V: VectorSpace<Field = f64> + Clone,
//     S: OdeSystem<Space = V>,
//     V: VectorSpace<Field = f64> + Clone,
// {
//     pub fn new() -> Self {
//         Self {
//             v: PhantomData,
//             s: PhantomData,
//         }
//     }

//     pub fn step(
//         &self,
//         f: S,
//         state: &TimeState<V>,
//         dt: f64, // <F::Space as VectorSpace>::Field,
//     ) -> TimeState<V> {
//         TimeState {
//             t: state.t + dt,
//             y: state.y.clone() + f(state) * dt,
//         }
//     }
// }

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
    V: VectorSpace + Clone,
    V::Field: Clone,
{
    fn solve_step(&self, f: &S, state: &TimeState<V>, dt: V::Field) -> TimeState<V> {
        self.step(f, state, dt)
    }
}

// impl OdeStepSolver for EulerSolver {
//     fn solve_step<S>(
//         &self,
//         f: &S,
//         state: &TimeState<S::Space>,
//         dt: <S::Space as VectorSpace>::Field,
//     ) -> TimeState<S::Space>
//     where
//         S: OdeSystem,
//         S::Space: Clone,
//         <S::Space as VectorSpace>::Field: Clone,
//     {
//         self.step(f, state, dt)
//     }
// }

// impl<V> OdeStepSolver for EulerSolver<V>
// where
//     V: VectorSpace<Field = f64> + Clone,
// {
//     fn solve_step<F>(
//         &self,
//         f: F,
//         state: &TimeState<F::Space>,
//         dt: <F::Space as VectorSpace>::Field,
//     ) -> TimeState<F::Space>
//     where
//         F: OdeSystem<Space = V>,
//     {
//         self.step(f, state, dt)
//     }
// }

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
            let ys = solver.integrate(&f, initial_state.clone(), t_end, n);

            let h = t_end / n as f64;
            let upper_bound = 20.0 * h;

            dbg!(k, n, h, &ys);

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
