import numpy as np

# inertial_coords = np.array([])
# rotating_coords = np.array([])

def inertial_to_rotating_frame(inertial_coords, rotation_angle):

    rotation_matrix = np.array([[np.cos(rotation_angle), -np.sin(rotation_angle), 0],
                                [np.sin(rotation_angle), np.cos(rotation_angle), 0],
                                [0, 0, 1]])


    rotating_coords = np.dot(rotation_matrix, inertial_coords)

    return rotating_coords


def rotating_to_inertial_frame(rotating_coords, rotation_angle):

    rotation_matrix = np.array([[np.cos(rotation_angle), np.sin(rotation_angle), 0],
                                [-np.sin(rotation_angle), np.cos(rotation_angle), 0],
                                [0, 0, 1]])


    inertial_coords = np.dot(rotation_matrix, rotating_coords)

    return inertial_coords
