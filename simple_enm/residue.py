"""
Contains the residue class.
"""
import numpy as np


class Residue:
    """Stores atributes of a protein residue"""

    def __init__(self, pdb_line):
        self.name = pdb_line[17:20].strip()
        self.chain = pdb_line[21:22].strip()
        self.residue_number = pdb_line[22:26].strip()
        self.coordinates = np.array(
            [
                float(pdb_line[30:38].strip()),
                float(pdb_line[38:46].strip()),
                float(pdb_line[46:54].strip()),
            ]
        )
        self.occupancy = float(pdb_line[54:60].strip())
        self.exp_b_factor = float(pdb_line[60:66].strip())

    def __str__(self):
        return f"{self.name} {self.chain} {self.residue_number} {self.coordinates} {self.occupancy} {self.exp_b_factor}"


if __name__ == "__main__":
    LINE = "ATOM      1  N   CYS A   1       5.142 -15.550  79.040  1.00 39.01           N  "
    cys = Residue(LINE)
    print(
        cys.name,
        cys.chain,
        cys.residue_number,
        cys.coordinates,
        cys.occupancy,
        cys.exp_b_factor,
    )
