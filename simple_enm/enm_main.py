import sys
import os
import numpy as np
from matplotlib import pyplot as plt
import motion_overlap as overlap


np.set_printoptions(precision=3)  # numpy settings
np.set_printoptions(suppress=True)
print("")


C = 2.99792458 * 10**8  # Speed of light


def normalize(array):
    """
    Takes a 1D array as input and returns a normalized version of the array
    """
    norm = np.linalg.norm(array)
    if norm == 0:
        return array
    return array / norm


def parse_pdb(file):
    """
    parses the pdb file. It considers only lines starting with 'ATOM' and atom
    name 'CA'.

    returns 1D arrays R0 and temp_factors. R0 contains the x, y, z
    coordinates of the C-alpha atoms. temp_factors contains the temperature
    facotr of each C-alpha atom. It returns also 3 lists containing the resnames
    , resids, names and chids.
    """
    R0 = list()
    exp_b = list()
    resnames = list()
    resids = list()
    names = list()
    chids = list()
    for line in file:
        if line.startswith("ATOM") and line[13:17].strip() == "CA":
            atom_serial = int(line[6:11].strip())
            atom_name = line[13:17].strip()
            residue = line[17:20].strip()
            chain = line[21:22].strip()
            res_id = line[22:26].strip()
            x_coord = float(line[30:38].strip())
            y_coord = float(line[38:46].strip())
            z_coord = float(line[46:54].strip())
            occup = float(line[54:60].strip())
            b = float(line[60:66].strip())
            if occup == 1.00:
                atom_coords = [x_coord, y_coord, z_coord]
                R0.extend(atom_coords)
                exp_b.append(b)
                names.append(atom_name)
                chids.append(chain)
                resnames.append(residue)
                resids.append(res_id)
    R0 = np.array(R0)
    exp_b = np.array(exp_b)
    return R0, exp_b, resnames, resids, names, chids


def parse_pdb_and_chain(file, chain_0):
    """
    parses the pdb file. It considers only lines starting with 'ATOM' and atom
    name 'CA'.

    returns 1D arrays R0 and temp_factors. R0 contains the x, y, z
    coordinates of the C-alpha atoms. temp_factors contains the temperature
    facotr of each C-alpha atom. It returns also 3 lists containing the resnames
    , resids, names and chids.
    """
    R0 = list()
    exp_b = list()
    resnames = list()
    resids = list()
    names = list()
    chids = list()
    for line in file:
        if (
            line.startswith("ATOM")
            and line[13:17].strip() == "CA"
            and line[21:22] == chain_0
        ):
            atom_serial = int(line[6:11].strip())
            atom_name = line[13:17].strip()
            residue = line[17:20].strip()
            chain = line[21:22].strip()
            res_id = line[22:26].strip()
            x_coord = float(line[30:38].strip())
            y_coord = float(line[38:46].strip())
            z_coord = float(line[46:54].strip())
            occup = float(line[54:60].strip())
            b = float(line[60:66].strip())
            if occup == 1.00:
                atom_coords = [x_coord, y_coord, z_coord]
                R0.extend(atom_coords)
                exp_b.append(b)
                names.append(atom_name)
                chids.append(chain)
                resnames.append(residue)
                resids.append(res_id)
    R0 = np.array(R0)
    exp_b = np.array(exp_b)
    return R0, exp_b, resnames, resids, names, chids


def build_distance_matrix(R0):
    """
    computes a NxN matrix with distances of different atoms. Takes as input
    an 1D array with 3N coordinates.

    returns a 2D NxN array
    """
    if R0.size % 3 != 0:
        print("The number of coordinates in R0 should be a multiple of 3!")
        return None
    n_atoms = int(R0.size / 3)
    distance_matrix = np.zeros((n_atoms, n_atoms), dtype=float)
    for i in range(n_atoms - 1):
        coords_i = R0[i * 3 : (i * 3) + 3]
        for j in range(i + 1, n_atoms):
            coords_j = R0[j * 3 : (j * 3) + 3]
            distance = np.linalg.norm(coords_i - coords_j)
            distance_matrix[i, j], distance_matrix[j, i] = distance, distance
    return distance_matrix


def build_sub_matrix(R0, distance_matrix, i, j, k=1):
    """
    Computes a 3x3 off diagonal sub-block of the mass weigted Hessian matrix
    that corresponds to atoms i and j. R0 is the array with 3N coordinates of
    all the atoms. distance_matrix contains all the distances between the atoms
    k is the force constant in the enm model.

    returns a 3x3 array or None if i is equal to j
    """
    if i == j:
        print("Division by zero")
        return None
    sub_matrix = np.empty((3, 3))
    coords_i = R0[i * 3 : (i * 3) + 3]
    coords_j = R0[j * 3 : (j * 3) + 3]
    coords_diff = coords_i - coords_j
    distance = distance_matrix[i, j]
    constant = (-k) * (1 / (distance**2))
    sub_matrix = constant * np.outer(coords_diff, coords_diff)
    return sub_matrix


def build_hessian(R0, Rc=13, k=1):
    """
    R0 is an array of x,y,z coordinates and Rc is an optional argument and it is
    the cutoff value in the enm model. k is the force constant in the enm model.
    The units Rc are angstrems and the units of k are kJ/A^2 mol.

    returns a 3Nx3N 2D ENM Hessian array.
    """
    n_atoms = int(R0.size / 3)
    hessian = np.zeros((R0.size, R0.size))
    distance_matrix = build_distance_matrix(R0)
    # Building the blocks (diagonals are different from off diagonal)
    for i in range(n_atoms):
        for j in range(i, n_atoms):
            if distance_matrix[i, j] > Rc:
                continue
            elif i != j:
                sub_matrix = build_sub_matrix(R0, distance_matrix, i, j, k=k)
                hessian[i * 3 : (i * 3) + 3, j * 3 : (j * 3) + 3] = sub_matrix
                hessian[j * 3 : (j * 3) + 3, i * 3 : (i * 3) + 3] = sub_matrix
            elif i == j:
                sub_matrix = np.zeros((3, 3))
                for atom in range(n_atoms):
                    if atom == i or distance_matrix[i, atom] > Rc:
                        continue
                    else:
                        sub_matrix += build_sub_matrix(
                            R0, distance_matrix, i, atom, k=k
                        )
                hessian[i * 3 : (i * 3) + 3, i * 3 : (i * 3) + 3] = (-1) * sub_matrix
    return hessian


def diagonalize(hessian):
    """
    computes the frequencies in cm-1 and eigenvectors of the hessian matrix.

    returns an 1D array of the frequencies in ascending order and an 2D array
    that contains the eigenvectors. The eigenvector corresponding to frequncy i
    is in the column i of the matrix.
    """
    e_val, e_vec = np.linalg.eigh(hessian)
    e_val_abs = abs(e_val) * 10**3  # Turning the eigenvalues into J kg-1 A-2
    freq = (np.sqrt(e_val_abs) / (2 * np.pi * C)) * 10**8  # 10^10 * 10^-2
    return e_val, freq, e_vec


def calculate_rmsf(freq, e_vec, T=1):
    """
    computes the rmsf in A-2 for all coordinates. freq is an array of
    frequencies, e_vec is a marix of eigenvectors. T is the temeprature.

    returns an array of rmsfs for each coordinate
    """
    N_coord = freq.size
    rmsf = np.zeros(N_coord)
    const = (8.3145 * T * 10 ** (20)) / (4 * np.pi**2)
    for i in range(N_coord):
        for k in range(6, N_coord):
            v = freq[k] * C * 10**2  # Converting cm-1 to s-1
            rmsf[i] = rmsf[i] + (e_vec[i, k] ** 2 / v**2)
    return rmsf * const


def calculate_temp_factors(rmsf):
    """
    calculates the temeprature factors for each atom and plots them together
    with the experimental values parsed form the pdb file.
    rmsf is the array of rmsf values for each coordinate and exp_b are the
    experimental temperature factor values for each atom.

    returns the array of calcualted temeprature_factors and plots them
    together with the expetimental values. Both values are normalized before
    plotting.
    """
    N = int(rmsf.size / 3)
    const = (8 * np.pi**2) / 3
    temp_factors = np.zeros(N)
    for i in range(N):
        temp_factors[i] = const * (np.sum(rmsf[i * 3 : i * 3 + 3]))
    temp_factors = normalize(temp_factors)
    return temp_factors


def create_b_plot(sub_out_folder, job_title, temp_factors, exp_b):
    """
    Creats and saves a plot of normalized experimental and calculated
    temperature factors B.
    """
    plot_file_name = job_title + "_b_plot.png"
    plot_path = os.path.join(sub_out_folder, plot_file_name)
    N = int(exp_b.size)
    exp_b = normalize(exp_b)
    x = np.arange(0, N)
    plt.plot(x, temp_factors, label="ENM theory")
    plt.plot(x, exp_b, label="Experiment")
    plt.legend(loc="best")
    plt.xlabel("Residue")
    plt.ylabel("Normalized temperature factor B")
    plt.savefig(plot_path)
    plt.close()


def create_correlation(sub_out_folder, job_title, temp_factors, exp_b):
    """
    Calculates the Pearson product-moment correlation coefficient between
    the calculated and experimental temeprature factors and saves the result
    inside a file.
    """
    out_file_name = job_title + "_corr" + ".txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    prs_corr_coeff = np.corrcoef(temp_factors, exp_b)[0, 1]
    with open(out_file_path, "w") as file:
        file.write(str(prs_corr_coeff))


def create_summary(
    sub_out_folder,
    pdb_model_id,
    pdb_target_id,
    n_residues,
    rmsd,
    overlap_arr,
    exp_b,
    temp_factors,
):
    """
    Creats a summary output file inside the sub_out_folder. It contains the
    following fields: pdb_mode, pdb_target, number of residues, RMSD, number of
    the mode with the biggest overlap, value of that overlap, cummulative
    overlap of the lowest 30 non zero vibrational modes, correlationa
    coefficient for the predicted B values.
    """
    out_file_name = "summary.txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    cum_overlap = np.sum(overlap_arr.copy()[0:36] ** 2)
    max_overlap = np.max(np.abs(overlap_arr)) ** 2
    max_overlap_mode = np.argmax(np.abs(overlap_arr))
    b_corr = np.corrcoef(temp_factors, exp_b)[0, 1]
    with open(out_file_path, "w") as file:
        file.write(f"{pdb_model_id},")
        file.write(f"{pdb_target_id},")
        file.write(f"{n_residues},")
        file.write(f"{rmsd:.5f},")
        file.write(f"{max_overlap_mode},")
        file.write(f"{max_overlap:.5f},")
        file.write(f"{cum_overlap:.5f},")
        file.write(f"{b_corr:.5f}")


def create_general_out(
    sub_out_folder, pdb_file_name, job_title, Rc, k, n_residues, e_val
):
    """
    Creats a general output file inside the sub_out_folder. If a file with
    the same name already exist it will get overwritten.
    """
    out_file_name = job_title + "_general.txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    pdb_id = pdb_file_name.removesuffix(".pdb")
    with open(out_file_path, "w") as file:
        file.write(f"PDB_id: {pdb_id}\n")
        file.write(f"Cutoff: {Rc}\n")
        file.write(f"Force_constant: {k}\n")
        file.write(f"Number of residues: {n_residues}\n")
        file.write("Small eigenvalues: ")
        np.savetxt(file, e_val[:20], fmt="%.3e", newline="  ")
        file.write("\n")


def create_nmd(
    names,
    resnames,
    chids,
    resids,
    exp_b,
    coordinates,
    e_vec,
    e_val,
    sub_out_folder,
    job_title,
):
    """
    The function creates a nmd file used for animating the normal modes in VMD
    The animation will be just for the first 10 modes. The exeprimental B
    factors are used for coloring the atoms that move the most.
    """
    out_file_name = job_title + "_vib" + ".nmd"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    title = pdb_file_name.replace(".pdb", ".anm")
    with open(out_file_path, "w") as file:
        file.write(f"nmwiz_load {out_file_name}\n")
        file.write(f"title {title}\n")
        file.write("names ")
        file.write(" ".join(names))
        file.write("\n")
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
        np.savetxt(file, coordinates, fmt="%.3f", newline=" ")
        file.write("\n")
        for i in range(6, 36):
            file.write(f"mode {i} {np.sqrt(1 / e_val[i])} ")
            np.savetxt(file, e_vec[:, i], fmt="%.3f", newline=" ")
            file.write("\n")


def create_eigenvalues_out(sub_out_folder, job_title, e_val):
    """
    Creats a file with all the eigenvalues inside the sub_out_folder.
    """
    out_file_name = job_title + "_eigenvalues.txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    with open(out_file_path, "w") as file:
        np.savetxt(file, e_val)


def create_eigenvectors_out(sub_out_folder, job_title, e_vec):
    """
    Creats a file with all the eigenvectors inside the sub_out_folder.
    """
    out_file_name = job_title + "_eigenvectors.txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    with open(out_file_path, "w") as file:
        np.savetxt(file, e_vec[:, :36], fmt="%9.5f")


def create_overlap_out(sub_out_folder, job_title, overlap_arr):
    """
    Creats a file with overlaps with the first 100 modes inside the
    sub_out_folder.
    """
    cum_overlap = np.sum(overlap_arr.copy()[0:36] ** 2)
    max_overlap = np.max(np.abs(overlap_arr)) ** 2
    out_file_name = job_title + "_indv_overlaps.txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    with open(out_file_path, "w") as file:
        file.write(f"CUM_OVERLAP(30): {cum_overlap}\n")
        file.write(f"MAX_OVERLAP: {max_overlap}\n")
        np.savetxt(file, overlap_arr, fmt="%9.5f")


def create_rmsd_out(sub_out_folder, job_title, rmsd, pdb_target_name):
    """
    Creats a file with the rmsd between two conformations.
    """
    out_file_name = job_title + "_rmsd.txt"
    out_file_path = os.path.join(sub_out_folder, out_file_name)
    with open(out_file_path, "w") as file:
        file.write(f"PDB target: {pdb_target_name}\n")
        file.write(f"RMSD: {rmsd}")


def create_cum_overlap_plot(sub_out_folder, job_title, overlap_arr):
    """
    Creats and saves a plot of cummulative overlaps of the first 30 normal
    modes.
    """
    plot_file_name = job_title + "_cum_overlap.png"
    plot_path = os.path.join(sub_out_folder, plot_file_name)
    N = int(overlap_arr.size)
    cum_overlap = np.cumsum(overlap_arr.copy() ** 2)
    x = np.arange(0, 36)
    plt.plot(x, cum_overlap[0:36])
    plt.title(job_title[0:4])
    plt.xlabel("Normal mode")
    plt.ylabel("Cummulative overlap")
    plt.savefig(plot_path)
    plt.close()


def create_overlap_plot(sub_out_folder, job_title, overlap_arr):
    """
    Creats and saves a plot of  overlaps of the first 100 modes.
    """
    plot_file_name = job_title + "_overlap.png"
    plot_path = os.path.join(sub_out_folder, plot_file_name)
    N = int(overlap_arr.size)
    x = np.arange(0, N)
    plt.plot(x, overlap_arr**2)
    plt.title(job_title[0:4])
    plt.xlabel("Normal mode")
    plt.ylabel("Individual overlap")
    plt.savefig(plot_path)
    plt.close()


# Parsing the input file (should be in the same folder as the python script)
input_file = sys.argv[1]
chain_flag = False  # Tracks if the chain is specified
with open(input_file, "r") as file:
    for line in file:
        split_line = line.split(maxsplit=1)
        if line.startswith("pdb_model"):
            pdb_model_id = split_line[1].strip()
            if len(split_line[1].strip()) == 4:
                pdb_file_name = f"{split_line[1].strip()}.pdb"
            elif len(split_line[1].strip()) == 5:
                chain_flag = True
                pdb_file_name = f"{split_line[1].strip()[:4]}.pdb"
                chain_0 = split_line[1].strip()[4].capitalize()
        if line.startswith("job_title"):
            job_title = split_line[1].strip()
        if line.startswith("pdb_folder"):
            pdb_folder = split_line[1].strip()
        if line.startswith("output_folder"):
            out_folder = split_line[1].strip()
        if line.startswith("cutoff"):
            Rc = float(split_line[1])
        if line.startswith("Force-constant"):
            k = float(split_line[1])
        if line.startswith("Temp-factors"):
            calc_b = split_line[1].strip()
        if line.startswith("overlap"):
            overlap_flag = split_line[1].strip()
        if line.startswith("pdb_target"):
            pdb_target_id = split_line[1].strip()
            if len(split_line[1].strip()) == 4:
                pdb_target_name = f"{split_line[1].strip()}.pdb"
            elif len(split_line[1].strip()) == 5:
                chain_flag = True
                pdb_target_name = f"{split_line[1].strip()[:4]}.pdb"
                chain_1 = split_line[1].strip()[4].capitalize()

print(chain_flag)
print(pdb_file_name)

# Creating a subdirectory inside the general output folder
sub_out_folder = os.path.join(out_folder, job_title)
try:
    os.mkdir(sub_out_folder)
except FileExistsError:
    pass

# Parsing the pdb file
file_path = os.path.join(pdb_folder, pdb_file_name)
with open(file_path, "r") as file:
    if chain_flag == True:
        R0, exp_b, resnames, resids, names, chids = parse_pdb_and_chain(file, chain_0)
    elif chain_flag == False:
        R0, exp_b, resnames, resids, names, chids = parse_pdb(file)
n_residues = len(resnames)

# Computing the eigenvalues and eigenvvectors of the Hessian
hessian = build_hessian(R0, Rc=Rc, k=k)
e_val, freq, e_vec = diagonalize(hessian)

# Computing the temperature factors, saving the plot, saving the corr_coeff
if calc_b == "True":
    rmsf = calculate_rmsf(freq, e_vec)
    temp_factors = calculate_temp_factors(rmsf)
    create_b_plot(sub_out_folder, job_title, temp_factors, exp_b)
    create_correlation(sub_out_folder, job_title, temp_factors, exp_b)

# Creating the output files
create_general_out(sub_out_folder, pdb_file_name, job_title, Rc, k, n_residues, e_val)
create_eigenvalues_out(sub_out_folder, job_title, e_val)
create_eigenvectors_out(sub_out_folder, job_title, e_vec)

# Creating the nmd file
create_nmd(
    names, resnames, chids, resids, exp_b, R0, e_vec, e_val, sub_out_folder, job_title
)

# Calculating the overlap
if overlap_flag == "True":
    pdb_model_path = os.path.join(pdb_folder, pdb_file_name)
    pdb_target_path = os.path.join(pdb_folder, pdb_target_name)
    if chain_flag == True:
        y, x = overlap.pdb_coords_and_chain(
            pdb_model_path, chain_0, pdb_target_path, chain_1
        )
    elif chain_flag == False:
        y, x = overlap.pdb_coords(pdb_model_path, pdb_target_path)
    nan_indices = overlap.get_nan_index(x)
    print(nan_indices)
    x_mod = np.delete(x, nan_indices)  # Removing the unmatched residues
    y_mod = np.delete(y, nan_indices)
    e_vec_mod = np.delete(e_vec, nan_indices, axis=0)
    x_mod = overlap.align(x_mod, y_mod)
    rmsd = overlap.calc_rmsd(x_mod - y_mod)
    overlap_arr = overlap.calc_all_overlaps(x_mod - y_mod, e_vec_mod)
    # Saving the overlap with the first 100 vibrational modes
    create_rmsd_out(sub_out_folder, job_title, rmsd, pdb_target_name)
    create_overlap_out(sub_out_folder, job_title, overlap_arr)
    create_overlap_plot(sub_out_folder, job_title, overlap_arr)
    create_cum_overlap_plot(sub_out_folder, job_title, overlap_arr)
    create_summary(
        sub_out_folder,
        pdb_model_id,
        pdb_target_id,
        n_residues,
        rmsd,
        overlap_arr,
        exp_b,
        temp_factors,
    )
