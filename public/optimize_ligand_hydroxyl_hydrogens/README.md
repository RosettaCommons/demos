# Optimize Ligand Hydroxyl Hydrogens

KEYWORDS: LIGANDS UTILITIES

This tutorial assumes a unix-style command line (or cygwin on windows).

## Build your ligand
In general the ligand atoms are marked by "HETATM", but check the ligand with pymol after you've grepped out the HETATM lines:
```
grep HETATM starting_inputs/cel5A_glucan.pdb > starting_inputs/cel5A_lig_noH.pdb
```
For this step you need to add on hydrogens. We use avogadro because its open source, but you can choose other software to place hydrogens:
http://avogadro.openmolecules.net
(this example is with 1.0.3 on mac)

Open in avogadro, choose Build--> add hydrogens, then save the molecule in MDL SDfile format cel5A_lig.mol

Alternatively, you can save the molecule in PDB format and convert the PDB file into a mol file using babel (http://openbabel.org/) as follows:
```
babel -ipdb starting_inputs/cel5A_lig.pdb -omol > starting_inputs/cel5A_lig.mol
```
Or the third alternative is to open the PDB file with pymol and save the molecule with a .mol extension to force pymol to save as a mol format. 

## Make a params file for the ligand
The next step is to create a params file for rosetta. Params file contains the internal coordinates of atoms, connectivity, charge of each atom, rosetta atom type. Most importantly for this demo it contains PROTON_CHI lines which specify the proton atoms that rosetta will simple around. 
-n specifies the name of the ligand in rosetta
```
python src/python/apps/public/molfile_to_params.py starting_inputs/cel5A_lig.mol -n cel
```
The output params file is called cel.params. 
The code will sample proton chi's at the explicit values stated in the params file, and then perform a minimization on those chi's.
You shouldn't have to change the params, but if you want to add sampling explicitly you can add more angles to sample. 
A sample proton chi is:
```
CHI 1  C1   C2   O1   H8
PROTON_CHI 1 SAMPLES 3 60 -60 180 EXTRA 0
```

## Setup files for rosetta
Pull out the protein without ligand:
```
grep ATOM starting_inputs/cel5A_glucan.pdb > rosetta_inputs/cel5A_input.pdb
```
Add back in the ligand pdb from molfile to params:
```
cat cel_0001.pdb >> rosetta_inputs/cel5A_input.pdb 
```

## Run rosetta enzdes
Now we can run the enzdes app in rosetta with minimal flags.
This optimizes proton chis on the ligand while also repacking sidechains:
```
path/to/EnzdesFixBB.[platform][compiler][mode] -s rosetta_inputs/cel5A_input.pdb -extra_res_fa cel.params -database path/to/minirosetta_database/ -out:file:o cel5A_score.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 10.0 -cut4 12.0 -minimize_ligand true
```
for example, using the provided inputs you can run: (where `$ROSETTA3`=path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/EnzdesFixBB.default.linuxgccrelease -s rosetta_inputs/cel5A_input.pdb -extra_res_fa rosetta_inputs/cel.params -out:file:o cel5A_score1.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 10.0 -cut4 12.0 -minimize_ligand true
```

Or, this optimizes proton chis on the ligand without repacking sidechains:
```
$> $ROSETTA3/bin/EnzdesFixBB.default.linuxgccrelease -s rosetta_inputs/cel5A_input.pdb -extra_res_fa rosetta_inputs/cel.params -out:file:o cel5A_score2.out -nstruct 1 -detect_design_interface -cut1 0.0 -cut2 0.0 -cut3 0 -cut4 0 -minimize_ligand true
```
Both runs should produce the PDB file cel5A_input__DE_1.pdb and the score file cel5A_score1.out or cel5A_score2.out, which can be placed into the output directory under different names, e.g. cel5A_output_nopack.pdb or cel5A_output_w_repack.pdb

Flag descriptions:
```
-extra_res_fa specifies the params file for new residues types (the glucan in this case)
-out:file:o specifies the file name for the enzdes-style score output. This contains extra information about the output design, like packing, interface energy, and many more
-nstruct 1 specifies one run and one output pdb; the packing is stochastic so for more sampling use a higher nstruct. 10-100 is recommended for most ligands.
-detect_design_interface tells rosetta to set up the designable and packable residues (in a packer task) based on distance from the ligand. Distances are calculated from every ligand heavy atom to the CA of amino acids.
[default values in brackets]
-cut1: CA less than cut1 is designed [6]
-cut2: CA between cut1 and cut2, with CA --> CB vector pointing towards ligand, is designed [8]
-cut3: CA less than cut3 is re-packed [10]
-cut4: CA between cut3 and cut4 with CA --> CB vector pointing towards ligand, is re-packed [12]
-minimize_ligand true  allow ligand torsions to minimize
```

## More Complete Energy Function / Sampling
If desired, use a more complete energy function and more sampling as in this example flags file
Recommended full flags for a more careful run:
```
enzdes_flags
```

## Full Enzyme Design (?)
This setup can also be used for full enzyme design
This run is very close to a full design of the active site. For a full design just change cut1 and cut2, e.g.
```
-cut1 6 -cut2 8
```
