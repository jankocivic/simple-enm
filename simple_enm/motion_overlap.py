"""
Module for calculating the overlap between normal modes and conformational
changes.
"""
import numpy as np


def normalize(array):
    """Normalizes an 1d array."""
    norm = np.linalg.norm(array)
    if norm == 0:
        return array
    return array / norm


def calc_center_of_mass(coordinate_matrix):
    """Calculates the center of mass of a system with coordinates of each
    residue stored in rows of a 2D numpy array.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return: vector of the center of mass
    """
    return coordinate_matrix.sum(axis=0) / coordinate_matrix.shape[0]


def shift_center_of_mass(coordinate_matrix):
    """Shifts the center of mass of a system to the origin.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return: a coordinate_matrix with the center of mass shifted to the origin
    """
    center_of_mass = calc_center_of_mass(coordinate_matrix)
    return coordinate_matrix - center_of_mass


def calc_corr_matrix(coord_matrix_0, coord_matrix_1):
    """Calculates the correlation matrix from two coordinate matrices. They
    should have the center of mass aligned.

    :param coord_matrix_0, coord_matrix_1: coordinate matrices whose rows are
        coordinates of each residue
    :return: the 3x3 correlation matrix
    """
    return coord_matrix_0.T @ coord_matrix_1


def calc_f_matrix(R):
    """Calculates a 4x4 F matrix from the correlation matrix R whose
    eigenvectors will correspod to the optimal rotational quaternion.

    :param R: 3x3 correlation matrix calculated by calc_corr_matrix
    :return f_matrix: explained above
    """
    f_matrix = [
        [
            R[0, 0] + R[1, 1] + R[2, 2],
            R[1, 2] - R[2, 1],
            R[2, 0] - R[0, 2],
            R[0, 1] - R[1, 0],
        ],
        [
            R[1, 2] - R[2, 1],
            R[0, 0] - R[1, 1] - R[2, 2],
            R[0, 1] + R[1, 0],
            R[0, 2] + R[2, 0],
        ],
        [
            R[2, 0] - R[0, 2],
            R[0, 1] + R[1, 0],
            R[1, 1] - R[0, 0] - R[2, 2],
            R[1, 2] + R[2, 1],
        ],
        [
            R[0, 1] - R[1, 0],
            R[0, 2] + R[2, 0],
            R[1, 2] + R[2, 1],
            R[2, 2] - R[1, 1] - R[0, 0],
        ],
    ]
    return f_matrix


def calc_rot_matrix(f_matrix):
    """Calculates the optimal rotational matrix from the eigenvector of the F
    matrix corresponding to the biggest eigenvalue.

    :param f_matrix: matrix generated by calc_f_matrix
    :return rot_matrix: the optimal rotational matrix
    """
    e_val, e_vec = np.linalg.eigh(f_matrix)
    q = e_vec[:, -1]
    rot_matrix = [
        [
            (q[0] ** 2 + q[1] ** 2 - q[2] ** 2 - q[3] ** 2),
            2 * (q[1] * q[2] - q[0] * q[3]),
            2 * (q[1] * q[3] + q[0] * q[2]),
        ],
        [
            2 * (q[1] * q[2] + q[0] * q[3]),
            (q[0] ** 2 - q[1] ** 2 + q[2] ** 2 - q[3] ** 2),
            2 * (q[2] * q[3] - q[0] * q[1]),
        ],
        [
            2 * (q[1] * q[3] - q[0] * q[2]),
            2 * (q[2] * q[3] + q[0] * q[1]),
            (q[0] ** 2 - q[1] ** 2 - q[2] ** 2 + q[3] ** 2),
        ],
    ]
    return rot_matrix


def rotate(coordinate_matrix, rot_matrix):
    """Rotates the system described by the coordinate_matrix.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :param rot_matrix: a 3x3 rotation matrix
    :return: rotated coordinate_matrix
    """
    rotated_coordinate_matrix = rot_matrix @ coordinate_matrix.T
    return rotated_coordinate_matrix.T


def align_coordinates(coord_matrix_0, coord_matrix_1):
    """Aligns the coordinates of coord_matrix_0 to the coordinates of
    coord_matrix_1 by minimizing the RMSD.

    :param coord_matrix_0, coord_matrix_1: coordinate matrices whose rows are
        coordinates of each residue
    :return: coord_matrix_0 after alignment
    """
    center_of_mass_1 = calc_center_of_mass(coord_matrix_1)
    coord_matrix_0 = shift_center_of_mass(coord_matrix_0)
    coord_matrix_1 = shift_center_of_mass(coord_matrix_1)
    R = calc_corr_matrix(coord_matrix_0, coord_matrix_1)
    F = calc_f_matrix(R)
    U = calc_rot_matrix(F)
    coord_matrix_0 = rotate(coord_matrix_0, U) + center_of_mass_1
    return coord_matrix_0


def calc_coordinate_difference(coord_matrix_0, coord_matrix_1):
    """Calculates the difference between coord_matrix_0 and coord_matrix_1.

    :param coord_matrix_0, coord_matrix_1: coordinate matrices whose rows are
        coordinates of each residue
    :return: flattend 1D coordinates difference vector
    """
    return (coord_matrix_0 - coord_matrix_1).flatten()


def calc_rmsd(coordinate_difference):
    """Calculates the rmsd from the coordinate difference vector of the aligned
    structures.

    :param coordinate_difference: 1D coordinates difference vector
    """
    rmsd = np.linalg.norm(coordinate_difference) / np.sqrt(
        (np.size(coordinate_difference) / 3)
    )
    return rmsd


def calc_overlaps(target_coord, model_coord, e_vec, n_modes=100):
    """Calculates the overlap between normal modes and conformational changes.

    :param target_coord: coordinate matrix of the conformer
    :param model_coord: coordinate matrix for the model system for which the
        normal modes were calculated
    :param e_vec: array where columns are individual eigenvectors is ascending
        order calculated for the model system
    :param n_modes: number of normal modes for which the overlap is calculated
    :return (overlaps, rmsd): touple containing an 1D array containing all
        individual overlaps and rmsd for the confomrational change
    """
    target_coord = align_coordinates(target_coord, model_coord)
    coordinate_difference = calc_coordinate_difference(target_coord, model_coord)
    n_modes = e_vec.shape()[1] if e_vec.shape()[1] <= n_modes else n_modes
    rmsd = calc_rmsd(coordinate_difference)
    return normalize(coordinate_difference) @ e_vec[:, :n_modes], rmsd
