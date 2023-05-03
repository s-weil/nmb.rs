use crate::ode_solvers::InitialValueOde;

// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

/*
public static double[] SecondOrder(double y0, double start, double end, int N, Func<double, double, double> f)
{
    double dt = (end - start) / (N - 1);
    double k1;
    double k2;
    double t = start;
    double[] y = new double[N];
    y[0] = y0;
    for (int i = 1; i < N; i++)
    {
        k1 = f(t, y0);
        k2 = f(t + dt, y0 + k1 * dt);
        y[i] = y0 + dt * 0.5 * (k1 + k2);
        t += dt;
        y0 = y[i];
    }
    return y;
} */

pub fn runge_kutta_second_order<const D: usize, F>(ode: F, t: f64, n: usize) -> Vec<[f64; D]>
where
    F: InitialValueOde<D>,
{
    let dt = t / n as f64;

    // let mut k1 = 0.0
    todo!()
}
