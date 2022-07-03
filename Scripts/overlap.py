"""
Script for calculating the overlap between normal modes and conformational
changes.
"""
import os
import sys
import numpy as np
from simple_enm import (
    pdb_tools,
    normal_modes,
    force_constant,
    temperature_factors,
    motion_overlap,
)


# INPUTS
JOB_TITLE = sys.argv[1]  # Used for naming the output file folder
MODEL_CODE = sys.argv[1]  # The fifth letter in the pdb code indicates one chain
TARGET_CODE = sys.argv[2]
FORCE_CONSTANT = "SIMPLE"  # For now SIMPLE is the only option
CUTOFF = 9  # Cutoff value in angstrems if the SIMPLE force constant was selected


def main():
    # Reading the pdb files
    model_path = os.path.join(os.getcwd(), f"PDB\\{MODEL_CODE[:4]}.pdb")
    target_path = os.path.join(os.getcwd(), f"PDB\\{TARGET_CODE[:4]}.pdb")
    with open(model_path, "r") as file:
        model_residues = pdb_tools.parse_pdb_file(file)
    with open(target_path, "r") as file:
        target_residues = pdb_tools.parse_pdb_file(file)

    # Matching the residue lists of the model and target
    pdb_tools.remove_partial_occup(model_residues)
    pdb_tools.remove_partial_occup(target_residues)
    if len(MODEL_CODE) == 5 and len(TARGET_CODE) == 5:
        model_residues = pdb_tools.select_chains(
            model_residues, chain_list=[MODEL_CODE[4]]
        )
        target_residues = pdb_tools.select_chains(
            target_residues, chain_list=[TARGET_CODE[4]]
        )
        (
            model_residues,
            target_residues,
            num_removed_residues,
        ) = pdb_tools.remove_unmatched(model_residues, target_residues)
    elif len(MODEL_CODE) == 4 and len(TARGET_CODE) == 4:
        (
            model_residues,
            target_residues,
            num_removed_residues,
        ) = pdb_tools.remove_unmatched(model_residues, target_residues, chain=False)
    else:
        raise ValueError(
            "Both the model and target PDB codes should be of the same length (4 or 5)"
        )

    # Calculation of the normal modes for the model system
    model_coordinates = normal_modes.build_coordinate_matrix(model_residues)

    distance_matrix = normal_modes.build_distance_matrix(model_coordinates)
    k_matrix = force_constant.k_with_cutoff(distance_matrix, cutoff=CUTOFF)

    hessian = normal_modes.build_hessian(model_coordinates, k_matrix, distance_matrix)
    e_val, e_vec = normal_modes.diagonalize_hessian(hessian)

    # Calculation of the temperature factors
    B_matrix = temperature_factors.calculate_temperature_factors(
        model_residues, e_val, e_vec
    )

    # Calculation of the overlaps
    target_coordinates = normal_modes.build_coordinate_matrix(target_residues)
    overlaps, rmsd = motion_overlap.calc_overlaps(
        target_coordinates, model_coordinates, e_vec, n_modes=1182
    )
    print(rmsd)
    print(np.cumsum(overlaps ** 2)[35])


if __name__ == "__main__":
    main()
