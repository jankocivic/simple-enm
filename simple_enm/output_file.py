"""
Contains functions for creating differnet output files.
"""

import os
import numpy as np
import matplotlib.pyplot as plt


# TXT OUTPUT FILES
def eigenvectors_out(out_folder_path, e_vec):
    """Creates an eigenvector.txt file that contains the first 36 eigenvectors
    inside the output folder for the job.

    :param out_file_path: path to the job subdirectory inside the Output folder
    :param e_vec: array where columns are individual eigenvectors in ascending
        order calculated for the model system
    """
    out_file_path = os.path.join(out_folder_path, "eigenvectors.txt")
    with open(out_file_path, "w") as file:
        np.savetxt(file, e_vec[:, :36], fmt="%9.5f")


def eigenvalues_out(out_folder_path, e_val):
    """Creates an eigenvalue.txt file that contains the eigenvalues in ascending
    order inside the output folder for the job.

    :param out_file_path: path to the job subdirectory inside the Output folder
    :param e_vec: array of eigenvalues in ascending order
    """
    out_file_path = os.path.join(out_folder_path, "eigenvalues.txt")
    with open(out_file_path, "w") as file:
        np.savetxt(file, e_val)


def enm_out(
    out_folder_path,
    pdb_file,
    chains,
    number_of_residues,
    force_constant,
    cutoff,
    b_corr,
):
    """To be added
    """
    out_file_path = os.path.join(out_folder_path, "general.txt")
    with open(out_file_path, "w") as file:
        file.write("Normal modes calculation\n")
        file.write(f"PDB_file: {pdb_file}\n")
        file.write(f"Chains: {chains}\n")
        file.write(f"Number of residues: {number_of_residues}\n")
        file.write(f"Force_constant: {force_constant}\n")
        file.write(f"Cutoff: {cutoff}\n")
        file.write(f"Temperature factor correlation: {b_corr}\n")


def overlap_out(
    out_folder_path,
    model_code,
    target_code,
    number_of_residues,
    number_of_removed,
    rmsd,
    force_constant,
    cutoff,
    b_corr,
    overlaps,
    e_val,
):
    """To be added
    """
    out_file_path = os.path.join(out_folder_path, "general.txt")
    overlaps_squared = overlaps[:100] ** 2
    max_overlap = np.max(overlaps_squared)
    max_overlap_mode = np.argmax(overlaps_squared)
    cum_overlap = np.sum(overlaps_squared[:36])
    with open(out_file_path, "w") as file:
        file.write("Overlap calculation\n")
        file.write(f"Model: {model_code}\n")
        file.write(f"Target: {target_code}\n")
        file.write(f"Number of residues: {number_of_residues}\n")
        file.write(f"Removed residues: {number_of_removed}\n")
        file.write(f"RMSD: {rmsd}\n")
        file.write(f"Force_constant: {force_constant}\n")
        file.write(f"Cutoff: {cutoff}\n")
        file.write(f"Temperature factor correlation: {b_corr}\n")
        file.write("Overlaps: ")
        np.savetxt(file, overlaps_squared, fmt="%9.5f", newline="  ")
        file.write("\n")
        file.write(f"max_overlap: {max_overlap}({max_overlap_mode})\n")
        file.write(f"cummulative ovarlap(30): {cum_overlap}")

    out_file_path = os.path.join(out_folder_path, "summary.txt")
    with open(out_file_path, "w") as file:
        file.write(f"{model_code},")
        file.write(f"{target_code},")
        file.write(f"{number_of_residues},")
        file.write(f"{number_of_removed},")
        file.write(f"{rmsd:.5f},")
        file.write(f"{e_val[6]},")
        file.write(f"{max_overlap_mode},")
        file.write(f"{max_overlap:.5f},")
        file.write(f"{cum_overlap:.5f},")
        file.write(f"{b_corr:.5f}")


# PLOTS
def b_plot(out_folder_path, pdb_code, temperature_factor_matrix):
    """To be added
    """
    plot_path = os.path.join(out_folder_path, "temperature_factors.png")
    calculated_b = temperature_factor_matrix[0]
    experimental_b = temperature_factor_matrix[1]
    number_of_residues = int(experimental_b.size)
    x_values = np.arange(0, number_of_residues)
    plt.title(pdb_code)
    plt.plot(x_values, calculated_b, label="ENM theory")
    plt.plot(x_values, experimental_b, label="Experiment")
    plt.legend(loc="best")
    plt.xlabel("Residue")
    plt.ylabel("Normalized temperature factor B")
    plt.savefig(plot_path)
    plt.close()


def overlap_plot(out_folder_path, model_code, target_code, overlaps):
    """To be added
    """
    plot_path = os.path.join(out_folder_path, "overlap.png")
    number_of_residues = int(overlaps.size)
    x_values = np.arange(0, number_of_residues)
    plt.plot(x_values, overlaps ** 2)
    plt.title(f"{model_code} --> {target_code}")
    plt.xlabel("Normal mode")
    plt.ylabel("Overlap")
    plt.savefig(plot_path)
    plt.close()


def cum_overlap_plot(out_folder_path, model_code, target_code, overlaps):
    """To be added.
    """
    plot_path = os.path.join(out_folder_path, "cum_overlap.png")
    cum_overlap = np.cumsum(overlaps.copy() ** 2)
    x_values = np.arange(0, 36)
    plt.plot(x_values, cum_overlap[0:36])
    plt.title(f"{model_code} --> {target_code}")
    plt.xlabel("Normal mode")
    plt.ylabel("Cummulative overlap")
    plt.savefig(plot_path)
    plt.close()


# NMD file
def create_nmd(
    out_folder_path, residue_list, pdb_file, e_vec, e_val,
):
    """To be added.
    """
    out_file_path = os.path.join(out_folder_path, "normal_modes.nmd")
    title = pdb_file.replace(".pdb", ".anm")
    resnames = [residue.name for residue in residue_list]
    chids = [residue.chain for residue in residue_list]
    resids = [residue.residue_number for residue in residue_list]
    exp_b = np.array([residue.exp_b_factor for residue in residue_list])
    coordinate_list = [residue.coordinates for residue in residue_list]
    coordinate_array = np.array(coordinate_list).flatten()
    with open(out_file_path, "w") as file:
        file.write("nmwiz_load normal_modes.nmd\n")
        file.write(f"title {title}\n")
        file.write("resnames ")
        file.write(" ".join(resnames))
        file.write("\n")
        file.write("chids ")
        file.write(" ".join(chids))
        file.write("\n")
        file.write("resids ")
        file.write(" ".join(resids))
        file.write("\n")
        file.write("betas ")
        np.savetxt(file, exp_b, fmt="%.2f", newline=" ")
        file.write("\n")
        file.write("coordinates ")
        np.savetxt(file, coordinate_array, fmt="%.3f", newline=" ")
        file.write("\n")
        for i in range(6, 36):
            file.write(f"mode {i} {np.sqrt(1 / e_val[i])} ")
            np.savetxt(file, e_vec[:, i], fmt="%.3f", newline=" ")
            file.write("\n")
