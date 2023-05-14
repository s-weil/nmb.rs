use crate::ode_solvers::{OdeState1D, OdeStepSolver1D};

pub struct EulerSolver;

impl EulerSolver {
    pub fn step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64,
    {
        let nexte_state = OdeState1D {
            t: state.t + dt,
            y: state.y + dt * f(state),
        };
        nexte_state
    }
}

impl OdeStepSolver1D for EulerSolver {
    fn solve_step<F>(&self, f: F, state: &OdeState1D, dt: f64) -> OdeState1D
    where
        F: Fn(&OdeState1D) -> f64,
    {
        self.step(f, state, dt)
    }
}
