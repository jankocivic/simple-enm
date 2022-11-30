"""
Script for calculating the normal modes of a protein. The script should be in a
directory that contains a "PDB" and "Output" subdirectory. The pdb files should
be located in the PDB directory and the output files will be saved in the Output
directory. The script also needs access to the simple_enm module. The default
cutoff is set to 12 A and can be changed manually.

It is necessary to specify two arguments:
:first argument: job title that will be used to name the output folder
:second argument: the PDB code of the protein

Example: "python enm_modes.py test1 1ALB"
"""

import os
import sys
# Locating the simple_enm module
sys.path.append("..\\")
sys.path.append("..\\simple_enm")

from simple_enm import (
    pdb_tools,
    normal_modes,
    force_constant,
    temperature_factors,
    output_file,
)


# Inputs
JOB_TITLE = sys.argv[1]
PDB_FILE = f"{sys.argv[2]}.pdb"
CHAINS = []
FORCE_CONSTANT = "SIMPLE"  # For now SIMPLE is the only option
CUTOFF = 12  # Cutoff value in angstrems if the SIMPLE force constant was selected


def main():
    pdb_file_path = os.path.join(os.getcwd(), f"PDB\\{PDB_FILE}")
    with open(pdb_file_path, "r") as file:
        residue_list = pdb_tools.parse_pdb_file(file)

    residue_list = pdb_tools.select_chains(residue_list, chain_list=CHAINS)
    residue_list = pdb_tools.remove_partial_occup(residue_list)

    number_of_residues = len(residue_list)

    coordinate_matrix = normal_modes.build_coordinate_matrix(residue_list)
    distance_matrix = normal_modes.build_distance_matrix(coordinate_matrix)

    k_matrix = force_constant.k_with_cutoff(distance_matrix, cutoff=CUTOFF)

    hessian = normal_modes.build_hessian(coordinate_matrix, k_matrix, distance_matrix)
    e_val, e_vec = normal_modes.diagonalize_hessian(hessian)

    temperature_factor_matrix = temperature_factors.calculate_temperature_factors(
        residue_list, e_val, e_vec
    )
    temperature_factor_correlation = temperature_factors.calculate_correlation_coefficient(
        temperature_factor_matrix
    )

    # Creating a subdirectory inside the general output folder
    out_folder_path = os.path.join(os.getcwd(), f"Output\\{JOB_TITLE}")
    try:
        os.mkdir(out_folder_path)
    except FileExistsError:
        pass

    # Output files
    output_file.eigenvectors_out(out_folder_path, e_vec)
    output_file.eigenvalues_out(out_folder_path, e_val)
    output_file.enm_out(
        out_folder_path,
        PDB_FILE,
        CHAINS,
        number_of_residues,
        FORCE_CONSTANT,
        CUTOFF,
        temperature_factor_correlation,
    )
    output_file.b_plot(out_folder_path, sys.argv[2], temperature_factor_matrix)
    output_file.create_nmd(out_folder_path, residue_list, PDB_FILE, e_vec, e_val)


if __name__ == "__main__":
    main()
