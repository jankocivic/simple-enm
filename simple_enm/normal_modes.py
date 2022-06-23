"""
Contains functions for performing normal mode analysis
"""
import numpy as np
import numba
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


@numba.jit(nopython=True)
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


@numba.jit(nopython=True)
def build_hessian_subblock(coordinates_i, coordinates_j, k_ij, d_ij):
    if d_ij == 0:  # Only happens when residue i and j are the same
        hessian_subblock = np.zeros((3, 3))
    else:
        coord_difference_vector = coordinates_j - coordinates_i
        hessian_subblock = np.outer(
            coord_difference_vector, coord_difference_vector
        ) * (-k_ij / (d_ij ** 2))
    return hessian_subblock


@numba.jit(nopython=True)
def build_outer_supermatrix(coordinate_matrix, k_matrix, distance_matrix):
    """Takes the (Nx3) coordinate matrix and creates a new (NxN) supermatrix
    where each block [i * 3 : (i * 3 + 3), j * 3 : (j * 3 + 3)] is a (3x3) matrix created
    from the outer product of the coordinates difference vector by itself,
    where the coordinates difference vector is the difference between the coordinate
    vector of residue j and residue i.

    :param coordinate_matrix: matrix whose rows are coordinates of each residue
    :return outer_supermatrix: explained above
    """
    number_of_atoms = coordinate_matrix.shape[0]
    outer_supermatrix = np.zeros((number_of_atoms * 3, number_of_atoms * 3))
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
                outer_supermatrix[i * 3 : (i * 3 + 3), j * 3 : (j * 3 + 3)] = subblock
                outer_supermatrix[j * 3 : (j * 3 + 3), i * 3 : (i * 3 + 3)] = subblock
    # Calculating the diagonal blocks as sums of the offdiagonal blocks
    for i in range(number_of_atoms):
        supermatrix_row = outer_supermatrix[i * 3 : (i * 3 + 3)]
        outer_supermatrix[
            i * 3 : (i * 3 + 3), i * 3 : (i * 3 + 3)
        ] = -supermatrix_row.reshape(number_of_atoms, 3, 3).sum(axis=0)
    return outer_supermatrix


# Preliminary testing
if __name__ == "__main__":
    with open("1PKL.pdb", "r", encoding="utf-8") as file:
        residue_list_outer = pdb.parse_pdb_file(file)
        coordinate_matrix = build_coordinate_matrix(residue_list_outer)
        distance_matrix = build_distance_matrix(coordinate_matrix)
        k_matrix = force_constant.k_with_cutoff(distance_matrix)
        C = build_outer_supermatrix(coordinate_matrix, k_matrix, distance_matrix)
