"""
Set of tools for working with PDB files.
"""
from residue import Residue


def parse_pdb_file(pdb_file):
    """Parses the pdb_file. Only considers the CA atom of each residue.
    Only considers lines starting with "ATOM".

    :param pdb_file: the pdb file
    :return residue_list: list storing all residues
    """
    residue_list = []
    for line in pdb_file:
        if line.startswith("ATOM") and line[13:17].strip() == "CA":
            residue_list.append(Residue(line))
    return residue_list


def select_chains(residue_list, chain_list):
    """Selects only particular protein chains from the list of all residues.
    Modifies the list of residues and returns it. If the list of chains is empty,
    all chains are selected.

    :param residue_list: list of all residues created by parse_pdb_file
    :param chains: a list with the selected protein chains
    :return residue_list_mod: modified residue list
    """
    if chain_list:
        residue_list_mod = [
            residue for residue in residue_list if residue.chain in chain_list
        ]
        return residue_list_mod
    return residue_list


def remove_partial_occup(residue_list):
    """Removes residues with an partial occupancy from the list of residues.

    :param residue_list: list of all residues created by parse_pdb_file
    "return residue_list_mod: modified residue list
    """
    residue_list_mod = [residue for residue in residue_list if residue.occupancy == 1]
    return residue_list_mod


def remove_unmatched(residue_list_0, residue_list_1, chain=True):
    """Creates an interescentions between two residue lists in the correct order.
    Creates new lists that contain the matched residues. The chain flag
    determines if only one chain is present which affects the comparison.

    :param residue_list_0, residue_list_1: lists of all residues created by the
        parse_pdb_file function
    :param chain: bool that determines the comparison criterion
    :return (residue_list_0_mod, residue_list_1_mod, num_removed_residues): 
        modified lists of residues  that contain only those residues that are 
        present in both lists together with the number of removed residues from
        the first residue list
    """
    residue_list_0_mod = []
    residue_list_1_mod = []
    if chain:
        for residue_0 in residue_list_0:
            for residue_1 in residue_list_1:
                if residue_0.compare_within_chain(residue_1):
                    residue_list_0_mod.append(residue_0)
                    residue_list_1_mod.append(residue_1)
                    break
    else:
        for residue_0 in residue_list_0:
            for residue_1 in residue_list_1:
                if residue_0.compare_overall(residue_1):
                    residue_list_0_mod.append(residue_0)
                    residue_list_1_mod.append(residue_1)
                    break
    num_removed_residues = len(residue_list_0) - len(residue_list_0_mod)
    return (residue_list_0_mod, residue_list_1_mod, num_removed_residues)
