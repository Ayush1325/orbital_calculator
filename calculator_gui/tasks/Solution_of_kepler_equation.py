"""the module is providing the solution of the kepler equation using the newton method"""

from numpy import sqrt, cbrt, sin, cos, arcsin


def solve_cubic(a, c, d):
    """ solve the cubic"""
    p = c/a
    q = d/a
    k = sqrt(q**2/4 + p**3/27)
    return cbrt(-q/2 - k) + cbrt(-q/2 + k)


def fun(e, m):
    """ solve for initial E"""
    n = sqrt(5 + sqrt(16 + 9/e))
    a = n*(e*(n**2 - 1)+1)/6
    c = n*(1-e)
    d = -m
    s = solve_cubic(a, c, d)
    return n*arcsin(s)


def find(e, m, E):
    """solve for the keplar next iteration deficient"""
    return E - e*sin(E) - m


def solve_kepler_equation(e, m):
    """solve for the Eccentric anomoly"""
    E = fun(e, m)
    tolerance = 1e-10
    # Newton's method
    while abs(find(e, m, E)) > tolerance:
        E -= find(e, m, E)/(1 - e*cos(E))
    return E


if __name__ == '__main__':

    # keep the eccentricty e:(0 <= e < 1)
    # keep the Mean anomoly M:(0 <= M <= pi)
    print(solve_kepler_equation(0.123, 0.343))
