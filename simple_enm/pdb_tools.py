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
    Modifies the list of residues and returns it. All chains are selected
    if the list of chains is empty.

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
