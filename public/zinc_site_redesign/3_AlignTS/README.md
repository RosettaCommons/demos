Step 3.1:
KEYWORDS: METALS UTILITIES
Align the transition state model on to the zinc site in the "clean" protein PDB file. 

Input:
- Protein PDB file: 1A4L_clean_A.pdb
- Ligand PDB file: LG_0001.pdb

Usage: 

```bash
python $ROSETTA_TOOLS/zinc_site_redesign/align.py -f 1A4L_clean_A.pdb -l LG_0001.pdb
```

Output:
- A file named aligned_ligand.pdb which is the transformed co-ordinates of the ligand after alignment (protein is held fixed).

Notes:

> **MANUALLY SET THE ATOMNAMES USED TO ALIGN CORRECTLY IN THE ALIGN.PY FILE. ORDER DOES MATTER: ZINC SHOULD BE SPECIFIED FIRST.**

------------------------------------------------------------------------------------------------------------------------------------------
Step 3.2

Generate files for next steps in Rosetta

Input:

- Clean PDB file: 1A4L_clean_A.pdb
- Aligned ligand: aligned_ligand.pdb
- name of metal: ZN

Usage:

```bash
python $ROSETTA_TOOLS/zinc_site_redesign/generate_metal_cstfile.py -f 1A4L_clean_A.pdb -m ZN -a aligned_ligand.pdb 
```

Output:

- Protein-ligand match contraints in a file: constraint.cst
- PDB for Rosetta: rosetta_cst.pdb

Other options:

-n to custom name the PDB for Rosetta 
