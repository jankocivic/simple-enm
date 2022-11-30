Short tutorial for using the python ENM program written by Janko Civic as part
of the TCCM internship in the 2021/2022 academic year.

1. The directory containing the program should contain the following:
  - readme.txt
  - Output  (Directory where the output files will be saved)
  - PDB   (Directory where we store PDB files)
  - enm_main.py
  - motion_overlap.py
  - enm.in

2. The input file contains the following fields:
  Job_title  (The name of the job that will be used for the output files)
  PDB_file  (The name of the pdb file for bulding the ENM)
  PDB_folder  (The path to the directory containing the PDB files)
  Output_folder (The path to the directory where the outputfiles will be saved)
  Cutoff  (The value of the cutoff for building the ENM)
  Force-constant  (Value of the force constant, should be left at 1)
  Temp-factors  ("True" if you want to calculate the temp-factors)
  Overlap  ("True" if you want to calculate the overlap with an other PDB file)
  PDB_target  (The name of the PDB file for which the overlaps will be calculated)

  Example of an input file:
  Job_title: 2hka_test
  PDB_file:  2hka.pdb
  PDB_folder: \home\janko\ENM_janko\PDB
  Output_folder: \home\janko\ENM_janko\Output
  Cutoff: 9
  Force-constant: 1
  Temp-factors: True
  Overlap: True
  PDB_target: 1nep.pdb

3. The program requirers the following python modules: sys, os, numpy,
  matplotlib, math

3. To run the program cd to the directory containing the program files and
   run it by typing "python enm_main.py enm.in"

4. The following output files are generated
  - Output
    - Job_title
      - Job_title_general (General information about the calcualtion)
      - Job_title_corr.txt (Correlation coefficient between temperature factors)
      - Job_title_eigenvalues.txt (Normal mode eigenvalues)
      - Job_title_eigenvectors.txt (First 30 eigenvectors)
      - Job_title_eigenvectors.npy (Eigenvectors in .npy format)
      - Job_title_rmsd (RMSD between two protein conformations)
      - Job_title_indv_overlaps.txt (Overlap with the first 100 normal modes)

      - Job_title_b_plot.png (Plot of the calculated and experimental b factors)
      - Job_title_overlap.png (Plot of the individual overlaps)
      - Job_title_cum_overlap.png (Plot of the cumulative overlaps)
      
      - Job_title_vib.nmd (File for visualizing the normal modes)

5. To visualize the normal modes it is necessary to install VMD.
  Open VMD and go to Extensions->Analysis->Normal Mode Wizard
  Choose the option "Load NMD file" and open the generated .nmd file


Remarks:
1. PDB files:
  Should only contain one model (careful with NMR structures).
  All HETATM entries are ignored.
  Only Ca atoms with occupancy 1 are considered.
  When comparing two conformations the two PDB files need to have the same
  primary sequance and not have too many missing residues.
  The program can run for more than 10 minutes for proteins with couple of
  thousand residues.
2. Outputfiles:
  The .npy file containing all the eigenvectors can become quite large for big
  proteins.