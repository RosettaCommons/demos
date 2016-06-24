Here we search for additional catalytic residues to stabilize the TS using RosettaMatch. 
In this example, we are looking for a Q to make Hbond to the attacking nucleophile hydroxyl.
Match-style constraint block corresponding to this interaction is appended to the constraints file obtained in the previous step. 
KEYWORDS: METALS DESIGN
The following files are needed

1. A protein scaffold file (1A4LA_clean_r.pdb) which includes a description of the pre-existing metal-protein interactions in the REMARK field. The protein co-ordinates could be the minimized ones from the previous CstOpt step or the starting scaffold co-ordinates without minimization (as in this example).
2. Flags for the matching: general_matching.flags subs.flags scaf.flags
3. all.pos: file describing what positions in the scaffold are to be used for which constraint (default all positions for each constraint)         
4. constraints.cst: geometric constraints including the pre-existing metal site constraints, and the "new" desired interaction.
5. Parameters for the ligand (LG.params) and optionally ligand rotamers as a concatanated PDB file - in this example, all atoms except the phosphorous and oxygens bonded to it are made virtual in the LG.params file to avoid using the rotamer ensemble for speed reasons.

Usage:

```bash
<path_to_Rosetta_directory>/main/source/bin/match.default.linuxgccrelease @general_matching.flags @scaf.flags @subs.flags -linmem_ig 10 -in:file::s 1A4L_clean_A_r.pdb
```

In the above, ".default.linuxgccrelease" may need to be updated for your build, operating system, and compiler.

Output:
If any matches are found new pdb files will be return with the residue inserted, named  UM_*pdb

Notes: For more details on the setup see the [matcher documentation](https://www.rosettacommons.org/docs/latest/application_documentation/design/match).
