from numpy import sqrt, cbrt, pi, sin, cos, arcsin


def solve_cubic(a, c, d):
    assert a > 0 and c > 0
    p = c / a
    q = d / a
    k = sqrt(q**2 / 4 + p**3 / 27)
    return cbrt(-q / 2 - k) + cbrt(-q / 2 + k)


def fun(e, M):
    n = sqrt(5 + sqrt(16 + 9 / e))
    a = n * (e * (n**2 - 1) + 1) / 6
    c = n * (1 - e)
    d = -M
    s = solve_cubic(a, c, d)
    return n * arcsin(s)


def main(e: float, M: float) -> float:
    assert 0 <= e < 1
    assert 0 <= M <= pi

    def f(E):
        return E - e * sin(E) - M

    E = fun(e, M)
    tolerance = 1e-10

    # Newton's method
    while abs(f(E)) > tolerance:
        E -= f(E) / (1 - e * cos(E))
    return E
