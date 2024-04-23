"""
NAME - Aashish Ranjan Singh
ADMISSION NO - 20JE0009

NAME - Ankit Kumar Singh
ADMISSION NO - 20JE0147
"""


import numpy as np
import math

def planet_elements_and_sv(planet_id, year, month, day, hour, minute, second):
    """
    This function calculates the orbital elements and the state
    vector of a planet from the date (year, month, day)
    and universal time (hour, minute, second).
    """
    global mu
    deg = math.pi / 180

    def J0(year, month, day):
        """
        This function calculates the Julian date at 0 hr UT of the given date.
        """
        A = math.floor(year / 100)
        B = 2 - A + math.floor(A / 4)
        j0 = math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day + B - 1524.5
        return j0

    def julian_date(year, month, day, hour, minute, second):
        """
        This function calculates the Julian Date from the provided date and time.
        """
        j0 = J0(year, month, day)
        ut = (hour + minute / 60 + second / 3600) / 24
        jd = j0 + ut
        return jd

    def zero_to_360(x):
        """
        This function reduces the angle to within the range of 0-360 degrees.
        """
        if x >= 360:
            x = x - math.floor(x / 360) * 360
        elif x < 0:
            x = x - (math.floor(x / 360) - 1) * 360
        return x

    def kepler_E(e, M):
        """
        This function solves Kepler's equation for eccentric anomaly (E).
        """
        E0 = M
        E1 = E0 + (M - E0 + e * math.sin(E0)) / (1 - e * math.cos(E0))
        while abs(E1 - E0) > 1.0e-06:
            E0 = E1
            E1 = E0 + (M - E0 + e * math.sin(E0)) / (1 - e * math.cos(E0))
        return E1

    def sv_from_coe(coe):
        """
        This function calculates the state vector (r,v) from the classical orbital elements (coe).
        """
        h, e, RA, incl, w, TA, a, w_hat, L, M, E = coe

        # Semi-major axis (a) in km
        a = a * au

        # Eccentricity (e) correction
        if e < 1.0e-08:
            e = 1.0e-08

        # Eccentric anomaly (E) in radians
        E = E * deg

        # Calculate position and velocity vectors
        p = a * (1 - e ** 2)
        r = np.array([(math.cos(TA * deg) - e), (math.sin(TA * deg) * math.sqrt(1 - e ** 2)), 0]) * p / (1 + e * math.cos(TA * deg))
        v = np.array([-math.sin(TA * deg), (e + math.cos(TA * deg)) * math.sqrt(1 - e ** 2), 0]) * math.sqrt(mu / p)

        return r, v

    # Orbital elements (will be calculated)
    coe = np.zeros(11)
    # State vector (will be calculated)
    r = np.zeros(3)
    v = np.zeros(3)

    # Equations
    jd = julian_date(year, month, day, hour, minute, second)
    J2000_elements, cent_rates = planetary_elements(planet_id)
    t0 = (jd - 2451545) / 36525
    elements = J2000_elements + cent_rates * t0
    a, e, incl, RA, w_hat, L = elements

    # Equation 2.61
    h = math.sqrt(mu * a * (1 - e ** 2))

    # Reduce the angular elements to within the range 0-360 degrees
    RA = zero_to_360(RA)
    w_hat = zero_to_360(w_hat)
    L = zero_to_360(L)
    w = zero_to_360(w_hat - RA)
    M = zero_to_360(L - w_hat)

    # Algorithm 3.1 (for which M must be in radians)
    E = kepler_E(e, M * deg)

    # Equation 3.10 (converting the result to degrees)
    TA = zero_to_360(2 * math.atan(math.sqrt((1 + e) / (1 - e)) * math.tan(E / 2)) / deg)

    coe = [h, e, RA, incl, w, TA, a, w_hat, L, M, E / deg]

    # Algorithm 4.2 (for which all angles must be in radians)
    r, v = sv_from_coe(coe)

    return coe, r, v, jd


def planetary_elements(planet_id):
    """
    This function provides the J2000 orbital elements and their rates for all planets.
    """
    au = 149597871
    deg = math.pi / 180

    J2000_elements = np.array([
        [0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084],
        [0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973],
        [1.00000011, 0.01671022, -0.00005, -11.26064, 102.94719, 100.46435],
        [1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332],
        [5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438],
        [9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432],
        [19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218],
        [30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003],
        [39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881]
    ])

    cent_rates = np.array([
        [0.00000066, 0.00002527, -23.51, -446.30, 573.57, 538101628.29],
        [0.00000092, -0.00004938, -2.86, -996.89, -108.80, 210664136.06],
        [-0.00000005, -0.00003804, -46.94, -18228.25, 1198.28, 129597740.63],
        [-0.00007221, 0.00011902, -25.47, -1020.19, 1560.78, 68905103.78],
        [0.00060737, -0.00012880, -4.15, 1217.17, 839.93, 10925078.35],
        [-0.00301530, -0.00036762, 6.11, -1591.05, -1948.89, 4401052.95],
        [0.00152025, 10.00019150, -2.09, -1681.4, 1312.56, 1542547.79],
        [-0.00125196, 0.00002514, -3.64, -151.25, -844.43, 786449.21],
        [-0.00076912, 0.00006465, 11.07, -37.33, -132.25, 522747.90]
    ])

    J2000_coe = J2000_elements[planet_id, :]
    rates = cent_rates[planet_id, :]

    # Convert from AU to km
    J2000_coe[0] *= au
    rates[0] *= au

    # Convert from arcseconds to fractions of a degree
    rates[2:6] = rates[2:6] / 3600

    return J2000_coe, rates


def zero_to_360(x):
    """
    This function reduces the angle to within the range of 0-360 degrees.
    """
    if x >= 360:
        x = x - math.floor(x / 360) * 360
    elif x < 0:
        x = x - (math.floor(x / 360) - 1) * 360
    return x


def month_planet_names(month, planet_id):
    """
    This function converts the month and planet ID numbers into names for output.
    """
    month_names = ["", "January", "February", "March", "April", "May", "June", "July",
                   "August", "September", "October", "November", "December"]
    planet_names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

    month_name = month_names[month]
    planet_name = planet_names[planet_id]

    return month_name, planet_name


def main(planet_id, year, month, day, hour, minute, second):
    # Input data
    # planet_id = 2
    # year = 2003
    # month = 8
    # day = 27
    # hour = 12
    # minute = 0
    # second = 0

    # Universal gravitational parameter (mu) in km^3/s^2
    mu = 1.327124e11
    # Convert degree to radian
    deg = math.pi / 180
    # AU to km
    au = 149597871

    # Call the function
    coe, r, v, jd = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second)

    # Convert the planet_id and month numbers into names for output
    month_name, planet_name = month_planet_names(month, planet_id)
    
    return(planet_name,year,month_name,day,hour,minute,second,jd,coe[0],coe[1],coe[2],coe[3],coe[4],coe[5],coe[6],coe[7],coe[8],coe[9],coe[10],r,np.linalg.norm(r),v,np.linalg.norm(v))

    # Echo the input data and output the solution to the command window
    # print("---------------------------------------------------")
    # print("Example 8.7\n")
    # print("Input data:\n")
    # print(f"Planet: {planet_name}")
    # print(f"Year : {year}")
    # print(f"Month : {month_name}")
    # print(f"Day : {day}")
    # print(f"Hour : {hour}")
    # print(f"Minute: {minute}")
    # print(f"Second: {second}\n")
    # print(f"Julian day: {jd:.3f}\n")
    # print("Orbital elements:")
    # print(f"\nAngular momentum (km^2/s) = {coe[0]}")
    # print(f"Eccentricity = {coe[1]}")
    # print(f"Right ascension of the ascending node (deg) = {coe[2]}")
    # print(f"Inclination to the ecliptic (deg) = {coe[3]}")
    # print(f"Argument of perihelion (deg) = {coe[4]}")
    # print(f"True anomaly (deg)  = {coe[5]}")
    # print(f"Semimajor axis (km)  = {coe[6]}")
    # print(f"Longitude of perihelion (deg)  = {coe[7]}")
    # print(f"Mean longitude (deg)  = {coe[8]}")
    # print(f"Mean anomaly (deg) = {coe[9]}")
    # print(f"Eccentric anomaly (deg)  = {coe[10]}\n")
    # print("State vector:")
    # print(f"\nPosition vector (km) = {r}")
    # print(f"Magnitude  = {np.linalg.norm(r)}")
    # print(f"Velocity (km/s) = {v}")
    # print(f"Magnitude = {np.linalg.norm(v)}")
    # print("-----------------------------------------------")