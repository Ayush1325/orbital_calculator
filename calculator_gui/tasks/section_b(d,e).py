# SHIVALI GUPTA
# 20JE0908



# SECTION B
# ---------------------------D-------------------------------------
# •	Gibbs’ method of preliminary orbit determination


import numpy as np

def gibbs_method(r1: np.ndarray, r2: np.ndarray, r3: np.ndarray, mu: float = 398600.4418) -> tuple:
    """
    Calculates orbital elements using Gibbs' method of preliminary orbit determination.

    Parameters:
        r1 (np.ndarray): Position vector at time t1.
        r2 (np.ndarray): Position vector at time t2.
        r3 (np.ndarray): Position vector at time t3.
        mu (float): Standard gravitational parameter of the central body.

    Returns:
        tuple: Orbital elements (semi-major axis, eccentricity, inclination,
               longitude of ascending node, argument of periapsis, true anomaly).
    """
    # Calculate cross products
    N = np.cross(r1, r2)
    D = np.cross(r2, r3)
    S = np.cross(r1 + r2 + r3, r2)
    
    # Calculate magnitudes
    N_mag = np.linalg.norm(N)
    D_mag = np.linalg.norm(D)
    S_mag = np.linalg.norm(S)
    
    # Calculate vectors for orbit determination
    R = r2
    V = (S_mag / D_mag) * (D / np.linalg.norm(r2)) - N / N_mag
    
    # Calculate orbital elements
    r = np.linalg.norm(R)
    v = np.linalg.norm(V)
    vr = np.dot(R, V) / r
    h = np.cross(R, V)
    h_mag = np.linalg.norm(h)
    energy = 0.5 * v**2 - mu / r
    
    # Semi-major axis
    a = -mu / (2 * energy)
    
    # Eccentricity vector
    e_vector = (1 / mu) * ((v**2 - mu / r) * R - r * vr * V)
    e = np.linalg.norm(e_vector)
    
    # Inclination
    i = np.arccos(h[2] / h_mag)
    
    # Longitude of ascending node
    n_hat = np.cross([0, 0, 1], h)
    n = np.linalg.norm(n_hat)
    omega = np.arccos(n_hat[0] / n)
    if n_hat[1] < 0:
        omega = 2 * np.pi - omega
    
    # Argument of periapsis
    omega_hat = np.cross(h, e_vector) / np.linalg.norm(np.cross(h, e_vector))
    omega_p = np.arccos(np.dot(n_hat, omega_hat) / (n * np.linalg.norm(omega_hat)))
    if e_vector[2] < 0:
        omega_p = 2 * np.pi - omega_p
    
    # True anomaly
    nu = np.arccos(np.dot(e_vector, R) / (e * r))
    if vr < 0:
        nu = 2 * np.pi - nu
    
    return a, e, np.degrees(i), np.degrees(omega), np.degrees(omega_p), np.degrees(nu)


# Make sure to provide position vectors r1, r2 and r3 
#   in kilometers. The output will give you the orbital elements:
#  semi-major axis, eccentricity, inclination, longitude of ascending
#  node, argument of periapsis, and true anomaly. 



#KINSHUK VERMA
#20JE0487


#SECTION B

# --------------------------E-------------------------
# Solution of Lambert’s problem for  orbit determination

import numpy as np
from scipy.optimize import brentq

def lambert(r1: np.ndarray, r2: np.ndarray, dt: float, mu: float = 398600.4418) -> tuple:
    """
    Solves Lambert's problem for orbit determination.
    
    Parameters:
        r1 (numpy.ndarray): Initial position vector.
        r2 (numpy.ndarray): Final position vector.
        dt (float): Time of flight between r1 and r2.
        mu (float): Standard gravitational parameter of the central body.
        
    Returns:
        (v1, v2): Initial and final velocity vectors.
    """
    # Define tolerance for convergence
    tol = 1e-8
    
    # Define a function to solve for x
    def func(x):
        A = np.sqrt(r1.dot(r2)) * (1 - np.cos(np.sqrt(mu) * x)) / np.sqrt(mu)
        y = np.linalg.norm(r1) + np.linalg.norm(r2) + A * (x ** 2)
        if A > 0:
            return y - np.sqrt(mu) * dt
        else:
            return y + np.sqrt(-mu) * dt
    
    # Solve for x using Brent's method
    x0 = brentq(func, -4 * np.pi, 4 * np.pi, xtol=tol)
    
    # Compute z value
    A = np.sqrt(r1.dot(r2)) * (1 - np.cos(np.sqrt(mu) * x0)) / np.sqrt(mu)
    B = np.sqrt(mu) * (1 - np.cos(np.sqrt(mu) * x0)) / np.sqrt(mu)
    z = np.dot(r1, r2) / np.sqrt(np.linalg.norm(r1) * np.linalg.norm(r2)) * np.sqrt(mu) * x0
    
    # Compute velocities
    f = 1 - A / np.linalg.norm(r1)
    g = B / np.sqrt(mu)
    fdot = np.sqrt(mu) / (np.linalg.norm(r1) * np.linalg.norm(r2)) * (z * np.sin(np.sqrt(mu) * x0) - np.sqrt(mu) * x0) * (1 - np.cos(np.sqrt(mu) * x0)) + (1 - np.cos(np.sqrt(mu) * x0)) / np.linalg.norm(r2)
    v1 = (r2 - f * r1) / g
    v2 = (fdot * r2 - r1) / g
    
    return v1, v2



# This code defines a function lambert() that takes the initial and final
#  position vectors (r1 and r2), the time of flight (dt), and optionally the 
# gravitational parameter (mu) of the central body. It returns the initial 
# and final velocity vectors (v1 and v2).