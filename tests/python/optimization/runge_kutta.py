import numpy as np
import matplotlib.pyplot as plt


def heun(f, x0, t):
    """Heun's method to solve x' = f(x,t) with x(t[0]) = x0.
    Code copied from https://www.math-cs.gordon.edu/courses/mat342/python/diffeq.py
    """
    n = len(t)
    x = np.array([x0] * n)
    for i in range(n - 1):
        h = t[i+1] - t[i]
        k1 = h * f(x[i], t[i])
        k2 = h * f(x[i] + k1, t[i+1])
        x[i+1] = x[i] + (k1 + k2) / 2.0

    return x


def rk2a(f, x0, t):
    """Second-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.
    Code copied from https://www.math-cs.gordon.edu/courses/mat342/python/diffeq.py
    """

    n = len(t)
    x = np.array([x0] * n)
    for i in range(n - 1):
        h = t[i+1] - t[i]
        k1 = h * f(x[i], t[i]) / 2.0
        x[i+1] = x[i] + h * f(x[i] + k1, t[i] + h / 2.0)

    return x


def rk2b(f, x0, t):
    """Second-order Runge-Kutta method to solve x' = f(x,t) with x(t[0]) = x0.
    Code copied from https://www.math-cs.gordon.edu/courses/mat342/python/diffeq.py
    """

    n = len(t)
    x = np.array([x0] * n)
    for i in range(n - 1):
        h = t[i+1] - t[i]
        k1 = h * f(x[i], t[i])
        k2 = h * f(x[i] + k1, t[i+1])
        x[i+1] = x[i] + (k1 + k2) / 2.0

    return x


def f(x, t):
    return x * np.sin(t)


a, b = (0.0, 10.0)
x0 = -1.0

n = 32
t = np.linspace(a, b, n)
h = t[1] - t[0]
tol = 1e-6

# compute various numerical solutions
x_heun = heun(f, x0, t)
x_rk2 = rk2a(f, x0, t)
# x_rk4 = rk4(f, x0, t)
# x_pc4 = pc4(f, x0, t)
# t_rkf, x_rkf = rkf(f, a, b, x0, tol, 1.0, 0.01)  # unequally spaced t

# compute true solution values in equal spaced and unequally spaced cases
x = -np.exp(1.0 - np.cos(t))
# xrkf = -np.exp(1.0 - np.cos(t_rkf))

plt.plot(t, x, 'b-o', t, x_heun, 'g-o', t, x_rk2, 'r-o')
plt.xlabel('$t$')
plt.ylabel('$x$')
plt.title('Solutions of $dx/dt = x\,\sin t$, $x(0)=-1$ ($h = %4.2f$)' % h)
plt.legend(('Solution', 'Heun $O(h^2)$', 'Runge-Kutta $O(h^2)$'),
           loc='lower left')
plt.show()
