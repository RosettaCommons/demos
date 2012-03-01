Task: Perform local docking of a flexible peptide

Why you might want to do this: 1) Check for the presence of a docking funnel, which would indicate a true binding site
2) Refine results from a global docking simulation 3) Build alternative backbones for reengineering the interaction 4)
To examine possible dynamics of bound structures

Steps: 
1) Select Chains to use in FlexPepDock - B & D were selected because they define the key interaction surface we
wish to model.

To extract chains B and D from 1OU8.pdb, use the following awk command from the local_dock_ssrA_peptide_against_sspB directory: 
awk '{ if( $1 =="ATOM" && ($5 == "B" || $5 =="D") ) print }' starting_files/1OU8.pdb > input_files/1OU8_BD.clean.pdb

1OU8_BD.clean.pdb no longer contains crystollagraphic waters. This completes preparing the input for FlexPepDock, which cannot handle explicit water molecules.

2) Prepack the starting structure
In most cases the initial should be "prepacked" before running FlexPepDock.  Prepacking will avoid starting structures with high energy side chain conformations far from the protein-peptide interface that could artificially inflate calculated energies, adding additional noise to simulations. 

path/to/FlexPepDocking.[platform][compiler][mode] -database path/to/rosetta_database -flexpep_prepack -in:file:s input_files/1OU8_BD.clean.pdb -out:path:pdb input_files -out:path:score:output_files -out:file:scorefile prepack.sc -ex1 -ex2aro

Options:
database: Specify the path to the Rosetta database, required for any Rosetta simulation
in:file:s: path to the input PDB file that was prepared previously
flexpep_prepack: prepack the structure
out:path:pdb: directory to write the prepacked structure to
out:path:score: directory to scorefiles to
out:file:scorefile: name of the scorefile 
ex1: Use extra sub-rotamers for the chi1 angle (+/- 1 standard deviation from the mean chi1 angle for each rotamer)
ex2aro: Use extra sub-rotamers from the chi2 angle of aromatic sidechains (+/-1 standard deviation from the mean chi2 for each aromatic rotamer)

Output:
This produces a prepacked structure named input_files/1OU8_BD.clean_0001.pdb and a scorefile named output_files/prepack.sc

Note that the prepacked structure may appear to have clashes at the interface.  This is expected behavior since the prepack step optimizes each docking partner in the unbound state.

If the out:path:pdb or out:path:score options are not specified, the corresponding files will be written to the working directory.  If out:file:scorefile is unspecified, the scorefile will be named score.sc.  If you have a scorefile with the same name in the scorefile directory, the new scores will be appended to the file as opposed to that file being overwritten.  It is important to make sure that the prepacking score is separated from the refinement scores.

3a) Refinement: Minimal movement

Commandline:

path/to/FlexPepDocking.[platform][compiler][mode] -database path/to/rosetta_database -pep_refine -out:path:pdb output_files -out:path:score output_files -out:file:scorefile refinement.sc -s input_files/1OU8_BD.clean_0001.pdb -ex1 -ex2aro -nstruct 3

Note: The output from the prepack step is the input for refinement

New options:
pep_refine: Run the high-resolution refinement mode of FlexPepDock.
nstruct: The number of decoys to generate.  It is common to generate O(10^5) decoys for best results.

Output:
The above commandline produces 3 candidate models (decoys) that will be named output_files/1OU8_BD.clean_0001_000x.pdb and a scorefile named output_files/refinement.sc. For explanation of the score terms used by FlexPepDock, visit http://www.rosettacommons.org/manuals/archive/rosetta[version]_user_guide/app_flexpep_docking.html  

Remarks:
We have included a PyMOL session of 5 superimposed decoys from a test run.  You can see there is very little movement of the peptide.  This type of simulation may be particularly useful for generating starting structures for design.

3b) Refinement: Moderate movement

Commandline:

path/to/FlexPepDocking.[platform][compiler][mode] -database path/to/rosetta_database -lowres_preoptimize -pep_refine -out:path:pdb output_files -out:path:score output_files -out:file:scorefile refinement.sc -s input_files/1OU8_BD.clean_0001.pdb -ex1 -ex2aro -nstruct 3

New options:
lowres_preoptimize: performs a low-resolution (centroid mode) stage with larger perturbations before the all-atom refinement stage

Remarks:
As opposed to 3a, this may be a better approach for detecting a docking funnel.

3c) Refinement: So much movement it's hardly refinement anymore

To add additional flexibility, one can supply "expert" flags that change the step size of the perturbations of the peptide backbone.  These flags are enumerated and explained in the FlexPepDock documentation (http://www.rosettacommons.org/manuals/archive/rosetta[version]_user_guide/app_flexpep_docking.html)

Commandline:
path/to/FlexPepDocking.[platform][compiler][mode] -database path/to/rosetta_database -lowres_preoptimize -pep_refine -smove_angle_range 12
 -out:path:pdb output_files -out:path:score output_files -out:file:scorefile refinement.sc -s input_files/1OU8_BD.clean_0001.pdb -ex1 -ex2aro -nstruct 3

New options:
smove_angle_range: Defines the perturbations size of small/shear moves.

Additional tips:
1) To reduce extraneous text printed to the screen during the run, supply the following flag:
-mute core basic

2) If there is a known structure that should be used for RMSD calculations, supply the following flag:
-in:file:native path/to/native.clean.pdb
This will add additional columns to the scorefile for the RMSD information.