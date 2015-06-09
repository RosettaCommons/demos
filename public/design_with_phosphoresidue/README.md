PART 1: MAKE SURE THE PHOSPHORESIDUE IS READ PROPERLY

We create the parameter file for the phosphoresidue as follow:

STEP 1: EXTRACTION OF PHOSPHOTYROSINE RESIDUE FROM PDB

We use one model of the NMR ensembles of the 2lct.pdb by deleting other models and keep the first one (3lct.pdb refers to first model).

From the "3lct.pdb" file, extract the part containing the atoms of the phosphoresidues 342 and 346 to generate "PT1.pdb" and "PT2.pdb". You can use the sample command line

grep PTR 3lct.pdb | grep 342 | grep HETATM > PT1.pdb
grep PTR 3lct.pdb | grep 346 | grep HETATM > PT2.pdb

STEP 2: CREATE MOLFILE FOR PTR FROM PDB FILE

-Either use avogadro(free chemical software) or go to the adress:
www.molecular-networks.com/online_demos/convert_demo
to generate an "PT1.mdl" file from the "PT1.pdb" with the original geometry of the phosphoresidue.

STEP 3: CREATE THE *.PARAMS FILE 
-Copy the "PTR.mdl" file (molfile format) into folder

mini/src/python/apps/public/molfile_to_params.py PT1.mdl -n PT1
mini/src/python/apps/public/molfile_to_params.py PT2.mdl -n PT2

This script creates PT1.params and PT1_0001.pdb. You'll need to delete the older PTR residues in 3lct.pdb and replace them with the co-ordinates for generated new residues PT1_0001.pdb and PT2_0001.pdb (4lct.pdb refers to 3lct.pdb with new phoshoresidues)

Also, change the "HETATM" tags for phoshoresidues in the PDB file to "ATOM".

PART 2: DESIGN AROUND THE PHOSPHORESIDUE

-Using pyMol, identify the neighbouring residues, within the desired radius (5 A in our case, PRT.resfile)

-Write the resfile. Use the option ALLAA next to the sequence position to design the selected residue to all of posibble amino acids.

sample command line:

~/mini/bin/fixbb.default.macosgccrelease -s 4lct.pdb -database ~/minirosetta_database -extra_res_fa PT1.params PT2.params -resfile PTR.resfile

This will generate a single pdb file with designed residues around the phosphoresidues. The sample outputs have been copied to /outputs.

If you want more structures, then use the flag -nstruct 1000 or -nstruct 10,000 depending on your need.



