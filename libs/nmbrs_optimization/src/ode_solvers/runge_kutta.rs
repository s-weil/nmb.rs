use crate::ode_solvers::T;

// TODO: add D=1 case for scalar like odes;
// split in steps, maybe even have a trait for just solving one step

// see https://people.math.ethz.ch/~hiptmair/tmp/NUMODE11.pdf

// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

pub fn runge_kutta_second_order_1d<F>(f: F, y0: f64, t_end: f64, n: usize) -> Vec<(f64, f64)>
where
    F: Fn(T, f64) -> f64,
{
    let dt = t_end / n as f64;
    let mut ys = Vec::with_capacity(n + 1);

    let mut t = 0.0;
    let mut y_t = y0;
    ys.push((t, y_t));

    while t < t_end {
        // evaluate f at t, y_t
        let k1 = f(t, y_t);
        // approximate f(t + dt, y1) via Euler f(t + dt, y + dt * f(t , y))
        let k2 = f(t + dt, y_t + dt * k1);
        // approximate y1 via Euler but the slope at t replaced by the mean of the slopes at t and t+dt, that is with k1 and k2
        let y = y_t + 0.5 * dt * (k1 + k2);
        y_t = y;
        t += dt;
        ys.push((t, y));
    }

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
    t_end: f64,
    n: usize,
) -> Vec<[f64; D]>
where
    F: Fn(T, &[T; D]) -> [T; D],
{
    let dt = t_end / n as f64;
    let mut t = 0.0;

    let mut ys = Vec::with_capacity(n + 1);

    let mut y0 = y0;

    while t < t_end {
        let k1 = f(t, &y0);
        let k2 = f(t + dt, &add(&y0, &scalar(dt, &k1)));
        let y = add(&y0, &scalar(dt * 0.5, &(add(&k1, &k2))));
        y0 = y;
        t += dt;
        ys.push(y);
    }

    assert_eq!(ys.len(), n);
    ys
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

    #[test]
    fn runge_kutta_second_order_1d_convegence() {
        // initial value problem
        let f = |t: f64, y| y * t.sin();
        let y0 = -1.0;

        // solution
        let sol = |t: f64| -(1.0 - t.cos()).exp();

        let t_end = 2.0;

        for k in 5..10 {
            let n = 2_usize.pow(k);
            let ys = super::runge_kutta_second_order_1d(f, y0, t_end, n);

            let h = t_end / n as f64;
            let upper_bound = 5.0 * h.powi(2);

            for i in 0..n {
                let (t_i, y_i) = ys[i];
                let sol_i = sol(t_i);
                let err_i = (sol_i - y_i).abs();
                assert!(err_i <= upper_bound);
            }
        }
    }
}
