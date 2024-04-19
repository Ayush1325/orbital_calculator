# Problem Statement : calculation of the state vector (r, v) given the initial state vector (r 0 , v 0 ) and the time lapse 

import numpy as np
import math

def main(R0: [float], V0: [float], t: int) -> tuple[[float], [float]]:
    global mu 
    # magnitude 
    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)

    vr0 = np.dot(R0, V0) / r0

    alpha = (2/r0) - (v0**2)/mu 
    print(f"alpha : {round(alpha * (10**5), 4)} x 10^(-5) Km^(-1)")

    x = Kepler_U(t, r0, vr0, alpha)
    print(f"x : {round(x, 4)} km^(1/2)")

    f, g = f_and_g(x, t, r0, alpha)
    print(f"f : {round(f, 4)} && g : {round(g, 4)} s^(-1)")

    R = f*R0 + g*V0


    r = np.linalg.norm(R)
    print(f"r : {r} km")

    fdot, gdot = fdot_and_gdot(x, r, r0, alpha)
    print(f"fdot : {round(fdot, 4)} s^(-1) && gdot : {round(gdot, 4)}")

    V = fdot * R0 + gdot * V0
    print(V)

    return R, V

def stumpC(z):
    if z > 0:
        c = (1 - math.cos(math.sqrt(z))) / z
    elif z < 0:
        c = (math.cosh(math.sqrt(-z)) - 1) / (-z)
    else:
        c = 1/2
    return c
    

def stumpS(z):
    if z > 0:
        s = (math.sqrt(z) - math.sin(math.sqrt(z))) / (math.sqrt(z))**3
    elif z < 0:
        s = (math.sinh(math.sqrt(-z)) - math.sqrt(-z)) / (math.sqrt(-z))**3
    else:
        s = 1/6
    return s


def f_and_g(x, t, r0, alpha):
    global mu
    f = 1 - (x**2 / r0) * stumpC(alpha * (x**2))
    g = t - (1/math.sqrt(mu)) * (x**3)*stumpS(alpha * (x**2))
    return f, g


def fdot_and_gdot(x, r, r0, alpha):
    global mu 
    fdot = (math.sqrt(mu) / (r*r0)) * (alpha*(x**3)*stumpS(alpha * (x**2)) - x)
    gdot = 1 - (x**2 / r) * stumpC(alpha * (x**2))
    return fdot, gdot


def Kepler_U(t, r0, vr0, alpha):
    error = 10**(-8)
    nMax = 1000
    global mu
    x = math.sqrt(mu) * abs(alpha) * t
    ratio = 1
    itr = 0

    
    while abs(ratio) > error and itr < nMax:
        z = alpha * (x**2)
        f_xi = (r0*vr0/math.sqrt(mu))*(x**2)*stumpC(z) + (1-alpha*r0)*(x**3)*stumpS(z) + r0*x - math.sqrt(mu)*t
        fdot_xi = ((r0*vr0/math.sqrt(mu))*x*(1-alpha*(x**2)*stumpS(z))) + (1 - alpha*r0)*(x**2)*stumpC(z) + r0
        ratio = f_xi / fdot_xi
        x = x - ratio

    return x



# Input Data
# example data 
global mu 
# default value of mu = 398600
mu = 398600

r_0 = []
v_0 = []
t = int(input("Time : "))

for i in range(0, 3):
    val = float(input("Value of R0 : "))
    r_0.append(val)

for i in range(0, 3):
    val = float(input("Value of V0 : "))
    v_0.append(val)


R0 = np.asarray(r_0)
V0 = np.asarray(v_0)

# R0 = np.array([1600, 5310, 3800])    # in KM
# V0 = np.array([-7.350, 0.4600, 2.470])  # in Km/s
# t = 3200    # in second 

R, V = main(R0, V0, t)

print(f"Initial Position Vector : {R0[0], R0[1], R0[2]} km")
print(f"Initial Velocity Vector : {V0[0], V0[1], V0[2]} km/s")
print(f"Time : {t} second")
print(f"Final Position Vector : {R[0], R[1], R[2]} km")
print(f"Final Velocity Vector : {V[0], V[1], V[2]} km/s")

