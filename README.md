# Simple ENM
Short tutorial for using the python ENM program written by Janko Civic as part
of the TCCM internship in the 2021/2022 academic year.

Good resources to learn about the theory of Elastic Network Models:
1. Yves-Henri Sanejouand, "Elastic Network Models: Theoretical and Empirical Foundations", *Methods Mol. Biol.*, **2012**, 924, 601-616
2. Monique M. Tirion, "Large Amplitude Elastic Motions in Proteins from a Single-Parameter. Atomic Analysis", *Phys. Rev. Lett.*, **1996**, 77, 1905

## Setup
1. Clone the github repository using the following command `git clone https://github.com/jankocivic/simple-enm.git`
2. Navigate to the cloned directory
3. Install the necessary dependencies in your local virtual environment
  1. Using conda you can create the environemnt and install the dependencies in the same time with the following command `conda env create -f environment.yml` and activate the environment with `conda activate simple-enm`
  2. If you don't have conda, create an environment and install the dependencies with `pip install -r requirements.txt` 

## Performing the calculations
The program is able to calculate the normal modes of a single protein and compute overlaps between the normal modes and conformational changes. The necessary scripts are located in the Scripts folder. The necessary PDB files need to be saved in the PDB folder and the output files will be stored in the Output folder.

### Normal modes
1. Navigate to the Scripts directory
2. Save the necessary PDB files in the PDB folder
3. Run the calcualtion with the following command `python enm_modes.py TITLE PDB_ID`
  - TITLE is the name of the directory that will be created in the Output folder and where output files will be stored
  - PDB_ID is the 4 letter PDB identifier (read notes at the end about PDB files)
4. The following output files are generated:
  - eigenvalues.txt \- Eigenvalues in ascending order
  - eigenvectors.txt \- The first 36 eigenvectors
  - general.txt \- General information about the calcualtion
  - normal_modes.nmd \- Files for visualizing the normal modes
  - temperature_factors.png \- Plot of the normalized calcualted and experimental temeprature factors of residue

### Overalp of normal modes and conformational changes 
1. Navigate to the Scripts directory
2. Save the necessary PDB files in the PDB folder
3. Run the calcualtion with the following command `python overlap.py PDB_ID[CHAIN_ID] PDB_ID[CHAIN_ID]`
  - first argument: PDB code of the model protein for which normal modes are calucluated
  - second argument: PDB code of the conformer
  - The PDB codes consist of 4 characters, but it is possible to append a fifth
letter to specify an exact chain
  - Change of the cutoff or the name of the output folder needs to be done manually inside the script
4. The following output files are generated:
  - eigenvalues.txt \- Eigenvalues in ascending order
  - eigenvectors.txt \- The first 36 eigenvectors
  - general.txt \- General information about the calcualtion
  - normal_modes.nmd \- Files for visualizing the normal modes
  - temperature_factors.png \- Plot of the normalized calcualted and experimental temeprature factors of residue
  - cum_overlap.png \- Plot of the cummulative overlap of the first N lowest modes and the conformational change
  - overlap.png \- Individual overlaps of each normal mode with the conformational change
  - summary.txt \- csv file with the most important data

### Visualization
1. Install VMD
2. Open VMD and go to:
   1. Extensions
   2. Analysis
   3. Normal Mode Wizard
   4. Load NMD file
   5. Select the desired .nmd file
   
### Remarks
- PDB files should only contain one model (careful with NMR structures)
- All HETATM entries are ignored
- Only Ca atoms with occupancy 1 are considered
- When comparing two conformations the two PDB files need to have the same residue numbering
- Too many missing residues can result in a disconnetcted graph causing the computaiton to fail
- The program can run for more than 10 minutes for proteins with couple of thousand residues
