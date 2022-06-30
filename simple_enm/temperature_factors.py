"""
Module for calculating temeprature factors.
"""
import numpy as np


def normalize(array):
    """Normalizes an 1d array."""
    norm = np.linalg.norm(array)
    if norm == 0:
        return array
    return array / norm


def calculate_rmsf(e_val, e_vec):
    """Calculates rmsf values for each coordinate. It ommits all the scaling
    constants becuase the units of the eigenvalues are arbitrary.

    :param e_val: array of eigenvalues in ascending order
    :param e_vec: matrix of eigenvectors corresponding to each eigenvalue
    :return rmsf: array of rmsf values for each coordinate
    """
    non_zero_e_val = e_val[6:]
    non_zero_e_vec = e_vec[:, 6:] ** 2
    reciprocal_e_val_matrix = np.diag(1 / non_zero_e_val)
    rmsf = np.sum(non_zero_e_vec @ reciprocal_e_val_matrix, axis=1)
    return rmsf


def calculate_temperature_factors(residue_list, e_val, e_vec):
    """Calculates temperature factors for each residue by taking the average of
    the rmsf values of each coordinate corresponding to the residue. Both the
    experimental and caluclated temeprature factors are normalized.

    :param residue_list: list of residues for which temp factors are computed
    :param e_val: array of eigenvalues in ascending order
    :param e_vec: matrix of eigenvectors corresponding to each eigenvalue
    :return temperature_factor_matrix: 2d array with the first and second rows
            are the normalized experimental and calculated temperature factors
    """
    # Consider only the first 1000 non-zero normal modes
    if np.size(e_val) > 1006:
        e_val = e_val[:1006]
        e_vec = e_vec[:, :1006]

    rmsf = calculate_rmsf(e_val, e_vec)
    number_of_residues = len(residue_list)
    experimental_temperature_factors = normalize(
        np.array([residue.exp_b_factor for residue in residue_list])
    )
    computed_temperature_factors = normalize(
        rmsf.reshape(number_of_residues, 3).sum(axis=1) / 3
    )
    temeprature_factor_matrix = np.row_stack(
        (experimental_temperature_factors, computed_temperature_factors)
    )
    return temeprature_factor_matrix


def calculate_correlation_coefficient(temperature_factor_matrix):
    """Calculates the Pearson product-moment correlation coefficient between
    the normalized experimental and calculated temeperature factors.

    :param temperature_factor_matrix: 2d array with the first and second rows
            are the normalized experimental and calculated temperature factors
    """
    correlation_coefficient = np.corrcoef(temperature_factor_matrix)[0, 1]
    return correlation_coefficient
