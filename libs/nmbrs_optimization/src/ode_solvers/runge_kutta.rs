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

    ys
}

#[cfg(test)]
mod tests {

    /*
           /// <summary>
           /// Runge-Kutta second order method for first order ODE.
           /// </summary>
           [Test]
           public void RK2Test()
           {
               Func<double, double, double> ode = (t, y) => t + 2 * y * t;
               Func<double, double> sol = (t) => 0.5 * (Math.Exp(t * t) - 1);
               double ratio = double.NaN;
               double error = 0;
               double oldError = 0;
               for (int k = 0; k < 4; k++)
               {
                   double y0 = 0;
                   double[] y_t = RungeKutta.SecondOrder(y0, 0, 2, Convert.ToInt32(Math.Pow(2, k + 6)), ode);
                   error = Math.Abs(sol(2) - y_t.Last());
                   if (oldError != 0)
                       ratio = Math.Log(oldError / error, 2);
                   oldError = error;
                   //Console.WriteLine(string.Format("{0}, {1}", error, ratio));
               }
               Assert.AreEqual(2, ratio, 0.01);// Check error convergence order
           }
    */

    #[test]
    fn runge_kutta_second_order() {
        let ode_f = |t: f64, y: &[f64; 1]| super::add(&[t], &super::scalar(2.0 * t, y));
        let _sol = |t: f64| 0.5 * (t.exp() - 1.0);
        let _ratio = f64::NAN;
        let _error = 0.0;
        let _old_error = 0.0;
        for k in 0..4 {
            let y0 = [0.0];
            let y_t = super::runge_kutta_second_order(ode_f, y0, 2.0, 2.0_f64.powi(k + 6) as usize);
            dbg!(y_t);
            // error = (sol(2.0) - y_t.last()).abs();
            // if old_error != 0.0 {
            //     ratio = (old_error / error).log();
            // }
            // old_error = error;
            //Console.WriteLine(string.Format("{0}, {1}", error, ratio));
        }
        assert!(false);
    }
}
