This is an optional step in the protocol where we optimize the protein-metal interaction in a polyAla context (i.e. everything in the pocket except the TS model and metal-chelating residues are trimmed back to Ala, and the constraints energy is optimized by minimization. See enzdes documentation for more details.)

Input files:

1. rosetta_cst.pdb: contains remark lines to specify catalytic residues, co-ordinates of the TS model
2. constraint.cst: contains match-style constraints between metal site and TS model
3. con_rotamers.pdb.gz : contains rotamers for the TS model
4. LG.params: parameter file for the TS model, which includes charges and topology
   of TS model. Here the following lines is added to the end of the file
   PDB_ROTAMERS con_rotamers.pdb.gz
   which ensures that the rotamers of the ligand are used for the optimization
5. optcst.flags
   Parameters to rosetta for the optimization

Usage:

./min.sh 

Output files:
PDB file for minimized interface:
e.g. rosetta_cst__DE_1.pdb

Typically this minimization is repeated multiple times to get slightly different conformations (number is controlled by the -nstruct option)

