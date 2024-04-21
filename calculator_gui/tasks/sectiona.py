import math


def kepler_hyperbolic(e: float, m: float, tol=1e-9, max_iter=1000) -> float | None:
    # Initial guess for eccentric anomaly (E)
    eccentric_anomaly = m / e

    # using Newton's method
    for _ in range(max_iter):
        f = e * math.sinh(eccentric_anomaly) - eccentric_anomaly - m
        f_prime = e * math.cosh(eccentric_anomaly) - 1
        g = f / f_prime
        eccentric_anomaly = eccentric_anomaly - g

        # Check for convergence
        if abs(g) < tol:
            return eccentric_anomaly

    return None


def kepler_equation(e: float, m: float, tol=1e-9, max_iter=1000) -> float | None:
    # Initial guess for eccentric anomaly (E)
    eccentric_anomaly = m / e

    # using Newton's method
    for _ in range(max_iter):
        f = eccentric_anomaly - e * math.sin(eccentric_anomaly) - m
        f_prime = 1 - e * math.cos(eccentric_anomaly)
        g = f / f_prime
        eccentric_anomaly = eccentric_anomaly - g

        # Check for convergence
        if abs(g) < tol:
            return eccentric_anomaly

    return None
