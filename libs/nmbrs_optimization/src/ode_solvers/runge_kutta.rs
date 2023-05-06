use crate::ode_solvers::T;

// TODO: add D=1 case fpr scalar like odes;
// split in steps, maybe even have a trait for just solving one step

// see https://people.math.ethz.ch/~hiptmair/tmp/NUMODE11.pdf

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

pub fn runge_kutta_second_order_1d<F>(f: F, y0: f64, t: f64, n: usize) -> Vec<f64>
where
    F: Fn(T, f64) -> f64,
{
    let dt = t / n as f64;
    let mut t = 0.0;

    let mut ys = Vec::with_capacity(n + 1);

    let mut y0 = y0;
    ys.push(y0);

    for _i in 0..n {
        // evaluate f at t, y0
        let k1 = f(t, y0);
        // evaluate f at t + dt, y0 + k1 * dt
        let k2 = f(t + dt, y0 + k1 * dt);
        let y = y0 + dt * 0.5 * (k1 + k2);
        y0 = y;
        ys.push(y);
        t += dt;
    }

    dbg!(ys.len(), t, y0);
    ys
}

fn scalar<const D: usize>(scalar: f64, vec: &[f64; D]) -> [f64; D] {
    let mut v = *vec;
    for i in 0..D {
        v[i] *= scalar;
    }
    v
}

fn add<const D: usize>(vec1: &[f64; D], vec2: &[f64; D]) -> [f64; D] {
    let mut v = *vec1;
    for i in 0..D {
        v[i] += vec2[i];
    }
    v
}

pub fn runge_kutta_second_order<const D: usize, F>(
    f: F,
    y0: [f64; D],
    t: f64,
    n: usize,
) -> Vec<[f64; D]>
where
    F: Fn(T, &[T; D]) -> [T; D],
{
    let dt = t / n as f64;
    let mut t = 0.0;

    let mut ys = Vec::with_capacity(n + 1);

    let mut y0 = y0;

    for _i in 0..=n {
        let k1 = f(t, &y0);
        let k2 = f(t + dt, &add(&y0, &scalar(dt, &k1)));
        let y = add(&y0, &scalar(dt * 0.5, &(add(&k1, &k2))));
        y0 = y;
        ys.push(y);
        t += dt;
    }

    assert_eq!(ys.len(), n);
    ys
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    #[test]
    fn runge_kutta_second_order_1d_convegence_order() {
        let a = 1.0;
        let b = 1.0;
        let f = |t, y| a * y + b;
        let y0 = 0.0;

        let sol = |t: f64| 0.5 * ((t * t).exp() - 1.0);

        let t = 2.0;
        let sol_2 = sol(t);

        let mut ratio: f64 = f64::NAN;
        let mut error = 0.0;
        let mut old_error = 0.0;

        for k in 5..10 {
            let n = 2_usize.pow(k);
            let ys = super::runge_kutta_second_order_1d(f, y0, t, n);

            for i in 0..n {
                let h = t / n as f64;
                let sol_i = sol(i as f64 * h);
                let err_i = (sol(i as f64 * h) - ys[i]).abs();
                let upper_bound = 100.0 * h.powi(2);
                dbg!(err_i, sol_i, ys[i], i, n, h, upper_bound);
                assert!(err_i <= upper_bound);
            }
        }
    }
    // #[test]
    // fn runge_kutta_second_order_1d_convegence_order() {
    //     let f = |t, y| t + 2.0 * y * t;
    //     let y0 = 0.0;

    //     let sol = |t: f64| 0.5 * ((t * t).exp() - 1.0);

    //     let t = 2.0;
    //     let sol_2 = sol(t);

    //     let mut ratio: f64 = f64::NAN;
    //     let mut error = 0.0;
    //     let mut old_error = 0.0;

    //     for k in 5..10 {
    //         let n = 2_usize.pow(k);
    //         let ys = super::runge_kutta_second_order_1d(f, y0, t, n);

    //         for i in 0..n {
    //             let h = t / n as f64;
    //             let err_i = (sol(i as f64 * h) - ys[i]).abs();
    //             let upper_bound = 40.0 * k as f64 * h.powi(2);
    //             dbg!(err_i, i, n, h, upper_bound);
    //             assert!(err_i <= upper_bound);
    //             // assert_abs_diff_eq!(err_i, h.powi(2), epsilon = 0.01);
    //         }
    //         // error = (sol_2 - y_t[n - 1]).abs();

    //         // if error != 0.0 {
    //         //     ratio = (old_error / error).log2();
    //         // }
    //         // old_error = error;
    //     }
    //     // assert_abs_diff_eq!(2.0, ratio, epsilon = 0.001);
    // }
}
