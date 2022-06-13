# overlap.py

import numpy as np


np.set_printoptions(precision=3)  # numpy settings
np.set_printoptions(suppress=True)


def normalize(array):
    """
    Takes a 1D array as input and returns a normalized version of the array
    """
    norm = np.linalg.norm(array)
    if norm == 0:
        return array
    return array / norm


def calc_center_of_mass(a):
    """
    Takes as input an numpy array of size 3N where N is the number of atoms and
    returns the position of the center of mass.
    """
    a = a.reshape(int(a.size / 3), 3)
    cm_coord = a.sum(axis=0) / a.shape[0]
    return cm_coord


def shift_center_of_mass(a):
    """
    Takes as input an numpy array of size 3N where N is the number of atoms and
    shifts the center of mass to the origin.
    """
    a = a.reshape(int(a.size / 3), 3)
    cm_coord = a.sum(axis=0) / a.shape[0]
    a = a - cm_coord
    return a.flatten()


def calc_corr_matrix(x, y):
    """
    Calculates the correlation matrix R taking as input 2 sets of coordinates
    of all the atoms. Both sets should have the same size and be amultiple of
    3. Returns a 3 by 3 correlation matrix.
    """
    x = x.reshape(3, int(x.size / 3), order="F")
    y = y.reshape(3, int(y.size / 3), order="F")
    R = x @ y.T
    return R


def calc_f_matrix(R):
    """
    Takes as input a 3 by 3 correlation matrix and returns a 4 by 4 F matrix
    whose eigenvectors will correspond to rotational quaternions.
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
    """
    Takes as input the F matrix and returns the 3 by 3 rotation matrix
    constructed from the eigenvector of F with the biggest eigenvalue.
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


def rotate(x, rot_matrix):
    """
    Takes as input an array x of 3N coordinates that we want to rotate and the
    rotation matrix. Outputs an array of rotated coordinates.
    """
    x = np.reshape(x, (3, int(x.size / 3)), order="F")
    x = rot_matrix @ x
    return x.flatten(order="F")


def align(x, y):
    """
    Takes as input two sets of 3N coordinates that correspond to the same atoms
    and alignes the coordinates x with y by minimizing RMSD using quaternions.
    Returns the vector of the aligned coordinates x.
    """
    cm_y = calc_center_of_mass(y)
    x = shift_center_of_mass(x)
    y = shift_center_of_mass(y)
    R = calc_corr_matrix(x, y)
    F = calc_f_matrix(R)
    U = calc_rot_matrix(F)
    x = rotate(x, U).reshape((-1, 3)) + cm_y
    return x.flatten()


def arr_coords_diff(x, y):
    """
    Takes as input two arrays of aligned 3N coordinates and returns the
    coordinate difference vector x-y.
    """
    return x - y


def pdb_coords(pdb_0, pdb_1):
    """
    Generates the coordinates of the two pdb files.
    """
    coords_0, coords_1 = list(), list()
    res_ids = list()  # Contains strings that help us locate the same residues
    with open(pdb_0, "r") as f0:
        for line in f0:
            if (
                line.startswith("ATOM")
                and line[13:17].strip() == "CA"
                and float(line[54:60].strip()) == 1
            ):
                res_ids.append(line[17:26].strip())
                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())
                coords_0.extend([x_coord, y_coord, z_coord])

        # Finding the equivalent reside in the second pdb file
        for r_id in res_ids:
            res_found = False
            with open(pdb_1, "r") as f1:
                for line in f1:
                    if (
                        line.startswith("ATOM")
                        and line[13:17].strip() == "CA"
                        and float(line[54:60].strip()) == 1
                        and line[17:26].strip() == r_id
                    ):
                        x_coord = float(line[30:38].strip())
                        y_coord = float(line[38:46].strip())
                        z_coord = float(line[46:54].strip())
                        res_found = True
                        break
            if res_found == True:
                coords_1.extend([x_coord, y_coord, z_coord])
            elif res_found == False:
                coords_1.extend([np.nan, np.nan, np.nan])
    coords_0 = np.array(coords_0)
    coords_1 = np.array(coords_1)
    return coords_0, coords_1


def pdb_coords_and_chain(pdb_0, chain_0, pdb_1, chain_1):
    """
    Generates the coordinates of the two pdb files if the chain is specified.
    """
    coords_0, coords_1 = list(), list()
    res_seq = list()  # Contains the residue sequence number
    with open(pdb_0, "r") as f0:
        for line in f0:
            if (
                line.startswith("ATOM")
                and line[13:17].strip() == "CA"
                and line[21:22].strip() == chain_0
                and float(line[54:60].strip()) == 1
            ):
                res_seq.append(line[22:26].strip())
                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())
                coords_0.extend([x_coord, y_coord, z_coord])

        # Finding the equivalent reside in the second pdb file
        for r_seq in res_seq:
            res_found = False
            with open(pdb_1, "r") as f1:
                for line in f1:
                    if (
                        line.startswith("ATOM")
                        and line[13:17].strip() == "CA"
                        and line[21:22].strip() == chain_1
                        and float(line[54:60].strip()) == 1
                        and line[22:26].strip() == r_seq
                    ):
                        x_coord = float(line[30:38].strip())
                        y_coord = float(line[38:46].strip())
                        z_coord = float(line[46:54].strip())
                        res_found = True
                        break
            if res_found == True:
                coords_1.extend([x_coord, y_coord, z_coord])
            elif res_found == False:
                coords_1.extend([np.nan, np.nan, np.nan])
    coords_0 = np.array(coords_0)
    coords_1 = np.array(coords_1)
    return coords_0, coords_1


def get_nan_index(a):
    """
    Takes as input an 1D array a and returns array of indices that contain nan.
    """
    return np.flatnonzero(np.isnan(a))


def calc_overlap(coords_diff, k_vec):
    """
    Calculates and returns the overlap of the motion with the eigenvector of a
    specific mode.
    """
    k_vec = normalize(k_vec)
    overlap = np.dot(coords_diff, k_vec) / (
        np.linalg.norm(k_vec) * np.linalg.norm(coords_diff)
    )
    return overlap


def calc_all_overlaps(coords_diff, e_vec, n_modes=100):
    """
    Calculates the overlap between the first n_modes non zero modes. Returns
    a numpy array with the individual overlaps.
    """
    overlaps = np.empty(n_modes)
    for i in range(0, n_modes):
        overlaps[i] = calc_overlap(coords_diff, e_vec[:, i])
    return overlaps


def calc_rmsd(coords_diff):
    """
    Calculates the rmsd from the coordinate difference vector of the aligned
    structures.
    """
    rmsd = np.linalg.norm(coords_diff) / np.sqrt((np.size(coords_diff) / 3))
    return rmsd
