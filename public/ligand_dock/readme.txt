Rosetta treats all biological objects as a string of residues. Each residue has a params file describing its coordinates. To dock a new ligand, the ligand has to have a params file that Rosetta understands, so the first step of any ligand task is ligand setup, follow this demo:
ligand_setup

After the ligand is ready, you can proceed, for example, to docking the ligand into a pocket in a pdb structure:
ligand_dock