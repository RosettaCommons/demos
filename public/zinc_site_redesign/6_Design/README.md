Next step is to optimize the interaction between TS model and protein
- design process - here new residues will be inserted in the protein
by Rosetta. 

Inputs:

1. A file from the matching step
2. Constraint file with the new constraint added
3. Rotamer file for the TS model
4. Parameter file for the TS model
5. Flag file for the design enzdes.flags

Usage:

./run.sh UM_1_H15H17H214D295Q58_1A4L_clean_A_r_1A4L_clean_A_1.pdb

Output:

UM_1_H15H17H214D295Q58_1A4L_clean_A_r_1A4L_clean_A_1__DE_1.pdb
enz_score.out

