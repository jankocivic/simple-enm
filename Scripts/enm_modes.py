"""
Script for calculating the normal modes of a protein. The script should be in a
directory that contains a "PDB" and "Output" subdirectory. The pdb files should
be located in the PDB directory and the output files will be saved in the Output
directory.
"""

# Inputs
JOB_TITLE = "1ALB_test"
PDB_FILE = "1ALB.pdb"
CHAINS = ["A"]
FORCE_CONSTANT = "SIMPLE"  # For now SIMPLE is the only option
CUTOFF = 12  # Cutoff value in angstrems if the SIMPLE force constant was selected


from simple_enm import pdb_tools, normal_modes, force_constant, temperature_factors
import os
import sys


def main():
    pdb_file_path = os.path.join(os.getcwd(), f"PDB\\{PDB_FILE}")
    with open(pdb_file_path, "r") as file:
        residue_list = pdb_tools.parse_pdb_file(file)
        print(residue_list[1])

    pdb_tools.select_chains(residue_list, chain_list=CHAINS)

    coordinate_matrix = normal_modes.build_coordinate_matrix(residue_list)
    distance_matrix = normal_modes.build_distance_matrix(coordinate_matrix)

    k_matrix = force_constant.k_with_cutoff(distance_matrix, cutoff=CUTOFF)

    hessian = normal_modes.build_hessian(coordinate_matrix, k_matrix, distance_matrix)
    e_val, e_vec = normal_modes.diagonalize_hessian(hessian)

    B_matrix = temperature_factors.calculate_temperature_factors(
        residue_list, e_val, e_vec
    )
    print(temperature_factors.calculate_correlation_coefficient(B_matrix))


if __name__ == "__main__":
    main()
