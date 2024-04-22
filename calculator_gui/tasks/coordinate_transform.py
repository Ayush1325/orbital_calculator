import numpy as np

"""
This module performs coordinate transformations between Geocentric equatorial and perifocal frames.
Name    : Tella Hemanth
Admn No : 20JE1027
"""

def compute_transformation_matrix(i, omega, Omega):
    """
    Compute the transformation matrix for converting between geocentric equatorial and perifocal frames.

    Args:
        i (float): Inclination angle in degrees.
        omega (float): Argument of periapsis angle in degrees.
        Omega (float): Longitude of ascending node angle in degrees.

    Returns:
        numpy.ndarray: Transformation matrix.
    """
    # Convert angles from degrees to radians
    i = np.radians(i)
    omega = np.radians(omega)
    Omega = np.radians(Omega)

    # Compute rotation matrices
    R1 = np.array([[1, 0, 0],
                   [0, np.cos(i), np.sin(i)],
                   [0, -np.sin(i), np.cos(i)]])
    
    R3_omega = np.array([[np.cos(omega), np.sin(omega), 0],
                         [-np.sin(omega), np.cos(omega), 0],
                         [0, 0, 1]])

    R3_Omega = np.array([[np.cos(Omega), np.sin(Omega), 0],
                         [-np.sin(Omega), np.cos(Omega), 0],
                         [0, 0, 1]])

    # Compute transformation matrix
    transformation_matrix = np.dot(R3_omega, np.dot(R1, R3_Omega))

    return transformation_matrix

def geocentric_to_perifocal(pos_geo, vel_geo, i, omega, Omega):
    """
    Convert geocentric equatorial coordinates to perifocal coordinates.

    Args:
        pos_geo (numpy.ndarray): Position vector in geocentric equatorial coordinates.
        vel_geo (numpy.ndarray): Velocity vector in geocentric equatorial coordinates.
        i (float): Inclination angle in degrees.
        omega (float): Argument of periapsis angle in degrees.
        Omega (float): Longitude of ascending node angle in degrees.

    Returns:
        tuple: Position and velocity vectors in perifocal coordinates.
    """
    # Compute transformation matrix
    transformation_matrix = compute_transformation_matrix(i, omega, Omega)

    # Apply transformation matrix
    pos_peri = np.dot(transformation_matrix, pos_geo)
    vel_peri = np.dot(transformation_matrix, vel_geo)

    return pos_peri, vel_peri

def perifocal_to_geocentric(pos_peri, vel_peri, i, omega, Omega):
    """
    Convert perifocal coordinates to geocentric equatorial coordinates.

    Args:
        pos_peri (numpy.ndarray): Position vector in perifocal coordinates.
        vel_peri (numpy.ndarray): Velocity vector in perifocal coordinates.
        i (float): Inclination angle in degrees.
        omega (float): Argument of periapsis angle in degrees.
        Omega (float): Longitude of ascending node angle in degrees.

    Returns:
        tuple: Position and velocity vectors in geocentric equatorial coordinates.
    """
    # Compute transformation matrix
    transformation_matrix = compute_transformation_matrix(i, omega, Omega)

    # Apply transpose of transformation matrix
    pos_geo = np.dot(transformation_matrix.T, pos_peri)
    vel_geo = np.dot(transformation_matrix.T, vel_peri)

    return pos_geo, vel_geo

"""
    Example Usage:
        Example position and velocity vectors in perifocal equatorial coordinates
        pos_peri = np.array([6285, 3628.6, 0])  --> Position vector (km)
        vel_peri = np.array([-2.4913, 11.290, 0])     --> Velocity vector (km/s)
        i = 30    --> Inclination (degrees)
        omega = 60  --> Argument of periapsis (degrees)
        Omega = 40  --> Longitude of ascending node (degrees)

        Convert from perifocal to geocentric equatorial
        pos_geo, vel_geo = perifocal_to_geocentric(pos_peri, vel_peri, i, omega, Omega)
        print("\nGeocentric equatorial coordinates:")
        print("Position:", pos_geo)
        print("Velocity:", vel_geo)
"""
