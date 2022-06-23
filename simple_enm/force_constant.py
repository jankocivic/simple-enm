"""
Contains functions for calculating different force constant matrices.
"""
import numpy as np


def k_with_cutoff(distance_matrix, cutoff=12):
    """Creates a matrix of force constants (k) between residues based on the distance matrix.
    The force constant equals to 1 for residues closer than the cutoff and 0 otherwise.

    :param distance_matrix: matrix of distances between the residues
    :param cutoff: the cutoff distance for the force constants in Angstrems
    :return k_matrix: matrix of force constants between residues
    """
    k_matrix = np.where(distance_matrix <= cutoff, 1, 0)
    return k_matrix


if __name__ == "__main__":
    distance_matrix_out = np.array([[0, 15, 5], [15, 0, 9], [5, 9, 0]])
    k_matrix_out = k_with_cutoff(distance_matrix_out)
    print(k_matrix_out)
