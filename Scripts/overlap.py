"""
Script for calculating the overlap between normal modes and conformational
changes. The script should be in a directory that contains a "PDB" and "Output"
subdirectory. The pdb files should be located in the PDB directory and the
output files will be saved in the Output directory. Modify the global JOB_TITLE
variable to change the name format of the created output folder.

It is necessary to specify two arguments:
:first argument: PDB code of the model protein for which normal modes are
    calucluated.
:second argument: PDB code of the conformer
The PDB codes consist of 4 characters, but it is possible to append a fifth
letter to specify an exact chain.

Example: "python overlap.py 3CHEA 2IUZB"
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
    motion_overlap,
    output_file,
)


# INPUTS
JOB_TITLE = f"{sys.argv[1]}_simple_15"  # Used for naming the output file folder
MODEL_CODE = sys.argv[1]  # The fifth letter in the pdb code indicates one chain
TARGET_CODE = sys.argv[2]
FORCE_CONSTANT = "SIMPLE"  # For now SIMPLE is the only option
CUTOFF = 15  # Cutoff value in angstrems if the SIMPLE force constant was selected


def main():
    print("\n")
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
    number_of_residues = len(model_residues)

    ########## TEMPORARY ###########
    print(f"{MODEL_CODE} --> {TARGET_CODE}")
    print(f"Number of removed residues: {num_removed_residues}")

    # Calculation of the normal modes for the model system
    model_coordinates = normal_modes.build_coordinate_matrix(model_residues)

    distance_matrix = normal_modes.build_distance_matrix(model_coordinates)
    k_matrix = force_constant.k_with_cutoff(distance_matrix, cutoff=CUTOFF)

    hessian = normal_modes.build_hessian(model_coordinates, k_matrix, distance_matrix)
    e_val, e_vec = normal_modes.diagonalize_hessian(hessian)

    ########### TEMPORARY #############
    print(f"Lowest eigenvalue: {e_val[6]}\n")

    # Calculation of the temperature factors
    temperature_factor_matrix = temperature_factors.calculate_temperature_factors(
        model_residues, e_val, e_vec
    )
    temperature_factors_correlation = temperature_factors.calculate_correlation_coefficient(
        temperature_factor_matrix
    )

    # Calculation of the overlaps
    target_coordinates = normal_modes.build_coordinate_matrix(target_residues)
    overlaps, rmsd = motion_overlap.calc_overlaps(
        target_coordinates, model_coordinates, e_vec, n_modes=100
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
    output_file.overlap_out(
        out_folder_path,
        MODEL_CODE,
        TARGET_CODE,
        number_of_residues,
        num_removed_residues,
        rmsd,
        FORCE_CONSTANT,
        CUTOFF,
        temperature_factors_correlation,
        overlaps,
        e_val,
    )
    output_file.b_plot(out_folder_path, MODEL_CODE, temperature_factor_matrix)
    output_file.overlap_plot(out_folder_path, MODEL_CODE, TARGET_CODE, overlaps)
    output_file.cum_overlap_plot(out_folder_path, MODEL_CODE, TARGET_CODE, overlaps)
    output_file.create_nmd(
        out_folder_path, model_residues, f"{MODEL_CODE[:4]}.pdb", e_vec, e_val
    )


if __name__ == "__main__":
    main()
