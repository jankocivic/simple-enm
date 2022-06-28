"""
Contains functions for performing normal mode analysis
"""
import numpy as np
import pdb_tools as pdb
import force_constant


def build_coordinate_matrix(residue_list):
    """Creates a matrix that contains the coordinates of each residue in rows.

    :param residue_list: list of all residues created by parse_pdb_file
    :return coordinate_matrix: matrix whose rows are coordinates of each residue
    """
    coordinate_list = [residue.coordinates for residue in residue_list]
    coordinate_matrix = np.array(coordinate_list)
    return coordinate_matrix


def build_distance_matrix(coordinate_matrix):
    """Creates a matrix of distances between the residues.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return distance_matrix: matrix of distances between the residues
    """
    coordinate_matrix_tile = np.tile(
        coordinate_matrix.T, (coordinate_matrix.shape[0], 1, 1)
    )
    distance_matrix = np.sqrt(
        ((coordinate_matrix_tile - coordinate_matrix_tile.T) ** 2).sum(axis=1)
    )
    return distance_matrix


def build_hessian_subblock(coordinates_i, coordinates_j, k_ij, d_ij):
    """Creates the off-diagonal 3x3 subblock of the Hessian matrix.

    :param coordinates_i: coordinates of residue i
    :param coordinates_j: coordinates of residue j
    :param k_ij: the ij element of the force constant matrix
    :param d_ij: the distance between reisidue i and j
    :return hessian_subblock: the 3x3 Hessian subblock between i and j
    """
    if d_ij == 0:  # Only happens when residue i and j are the same
        hessian_subblock = np.zeros((3, 3))
    else:
        coord_difference_vector = coordinates_j - coordinates_i
        hessian_subblock = np.outer(
            coord_difference_vector, coord_difference_vector
        ) * (-k_ij / (d_ij ** 2))
    return hessian_subblock


def build_hessian(coordinate_matrix, k_matrix, distance_matrix):
    """Creates the Hessian matrix. Each off-diagonal subblock is calculated by
    the build_hessian_subblock function. The digaonal subblocks are negative
    sums of all off-diagonal subblocks in the same row.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return hessian: explained above
    """
    number_of_atoms = coordinate_matrix.shape[0]
    hessian = np.zeros((number_of_atoms * 3, number_of_atoms * 3))
    # Calculating offdiagonal blocks
    for i in range(number_of_atoms - 1):
        for j in range(i + 1, number_of_atoms):
            if k_matrix[i, j] != 0:
                subblock = build_hessian_subblock(
                    coordinate_matrix[i],
                    coordinate_matrix[j],
                    k_matrix[i, j],
                    distance_matrix[i, j],
                )
                hessian[i * 3 : (i * 3 + 3), j * 3 : (j * 3 + 3)] = subblock
                hessian[j * 3 : (j * 3 + 3), i * 3 : (i * 3 + 3)] = subblock
    # Calculating the diagonal blocks as sums of the offdiagonal blocks
    for i in range(number_of_atoms):
        hessian_row = hessian[i * 3 : (i * 3 + 3)]
        hessian[i * 3 : (i * 3 + 3), i * 3 : (i * 3 + 3)] = -sum(
            np.hsplit(hessian_row, number_of_atoms)
        )
    return hessian


def diagonalize_hessian(hessian):
    """
    Finds eigenvalues and eigenvectors of the hessian matrix. Just
    a wrapper for the numpy.linalg.eigh function.

    :param hessian: the Hessian matrix
    :return: array of eigenvalues in ascending order, array of eigenvectors
             where column [i] is the eigenvector with eigenvalue [i]
    """
    return np.linalg.eigh(hessian)


# Preliminary testing
if __name__ == "__main__":
    with open("1ALB.pdb", "r", encoding="utf-8") as file:
        residue_list_outer = pdb.parse_pdb_file(file)
        coordinate_matrix = build_coordinate_matrix(residue_list_outer)
        distance_matrix = build_distance_matrix(coordinate_matrix)
        k_matrix = force_constant.k_with_cutoff(distance_matrix)
        C = build_hessian(coordinate_matrix, k_matrix, distance_matrix)
        e_val, e_vec = diagonalize_hessian(C)
        np.savetxt("1ALB_hessian", C, fmt="%9.5f")
        np.savetxt("1ALB_evec.txt", e_vec, fmt="%9.5f")
        np.savetxt("1ALB_eval.txt", e_val)
