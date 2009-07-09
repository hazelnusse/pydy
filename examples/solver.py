from numpy import array, empty

def rk4_step(x, y, h, derivs):
    """
    Does one Runge-Kutta 4th order step.

    Use rk4int() to solve the equations on the whole interval.
    """
    k1 = h * derivs(y, x)
    k2 = h * derivs(y + k1/2, x + h/2)
    k3 = h * derivs(y + k2/2, x + h/2)
    k4 = h * derivs(y + k3, x + h)
    return y + k1/6 + k2/3 + k3/3 + k4/6

def rk4int(derivs, y0, t, h=0.1):
    """
    Solves the set of linear ODEs using rk4.

    Example:
    >>> t = arange(0, 7.e-5, h)
    >>> sol = rk4int(derivs, t, (a0, b0, c0))
    >>> from pylab import plot, show, legend
    >>> plot(t, sol[:, 0], ".", label="y_1")
    >>> plot(t, sol[:, 1], ".", label="y_2")
    >>> plot(t, sol[:, 2], ".", label="y_3")
    """
    y = array(y0)

    result = empty((len(t), len(y)))
    for i, x in enumerate(t):
        y = rk4_step(x, y, h, derivs)
        result[i, :] = y
    return result

if __name__ == "__main__":
    # simple example from chemical kinetics:
    from numpy import arange
    from scipy.integrate import odeint

    def derivs(y, x):
        k1 = 9.e5
        k2 = 1.e5
        k3 = 3.e5
        return array([
            -(k1 + k3) * y[0]            ,
                 k1    * y[0] - k2 * y[1],
                 k3    * y[0] + k2 * y[1],
            ])


    a0 = 9.e10
    b0 = 0.
    c0 = 6.e9

    h = 1e-7
    t = arange(0, 7.e-5, h)

    sol = rk4int(derivs, (a0, b0, c0), t)
    sol_scipy = odeint(derivs, (a0, b0, c0), t)

    from pylab import plot, show, legend
    plot(t, sol[:, 0], ".", label="A rk4")
    plot(t, sol[:, 1], ".", label="B rk4")
    plot(t, sol[:, 2], ".", label="C rk4")
    plot(t, sol_scipy[:, 0], "x", label="A scipy")
    plot(t, sol_scipy[:, 1], "x", label="B scipy")
    plot(t, sol_scipy[:, 2], "x", label="C scipy")
    legend()
    show()
