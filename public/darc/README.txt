--------------------------------------------------------------------------------
   INSTRUCTIONS FOR RUNNING DARC PROTOCOL
--------------------------------------------------------------------------------

	The Docking Approach using Ray Casting (DARC) protocol is used to dock the 
small molecule in a pocket centered around residue 105 of the protein, B-cell 
lymphoma-extra large, 2YXJ.pdb. 

--------------------------------------------------------------------------------
   GENERATING INPUT FILES
--------------------------------------------------------------------------------

	Input files are included in the folder inputs, but the following explains 
how they were obtained or generated. 2YXJ.pdb was downloaded from The Protein Databank.
The small molecule used for docking was found in ZINC, the database of commercially 
available compounds at http://zinc.docking.org/. We download ZINC00057615 in SMILES 
format. This was converted to a pdb file with 100 conformers using OMEGA2 from open eye, 
http://www.eyesopen.com/ using the following command:

$ OpenEye/bin/omega2 -n ZINC00057615.smi -out molecules.pdb -maxconformers 100

	Next, the file with all the conformers is split into multiple files with one 
conformer per file with babel from python using the following command:

$ babel -pdb molecules.pdb -opdb molecule.pdb -m

	Next, the sdf files are converted to pdb and parameter files using molfile_to_params 
which is a python app in rosetta. This command is used:

$ babel -pdb molecule.pdb -opdb molecule.sdf
$ path/molfile_to_params.py -c -nKHR -pmol molecule.sdf

--------------------------------------------------------------------------------
   BUILD ROSETTA
--------------------------------------------------------------------------------

$ scons mode=release bin

To build with GPU enabled:
$ scons mode=release extras=opencl bin

--------------------------------------------------------------------------------
   Run DARC
--------------------------------------------------------------------------------

$ cd run_dir	

Copy 2YXJ.pdb, molecule.pdb, and molecule.params to run_dir. 
Run DARC with the following command:
$ path/DARC.opencl.linuxgccrelease -database ~/rosetta_database -input_protein_file 2YXJ.pdb 
  -input_ligand_file molecule.pdb -extra_res_fa molecule.params -gpu 0 -eggshell_triplet rays.txt

Run DARC with the following command using GPU:
$ path/DARC.opencl.linuxgccrelease -database ~/rosetta_database -input_protein_file 2YXJ.pdb 
  -input_ligand_file molecule.pdb -extra_res_fa molecule.params -gpu 1 -eggshell_triplet rays.txt

Output files:
Running DARC will result in this output file:
DARC_molecule_2YXJ_-1.pdb is the final docked pose

--------------------------------------------------------------------------------
   ANALYZING THE RESULTS
--------------------------------------------------------------------------------

The output should look similar to:

$ ~/rosetta_source/bin/DARC.opencl.linuxgccrelease -database ~/rosetta_database -input_protein_file 2YXJ.pdb 
  -input_ligand_file molecule.pdb -extra_res_fa molecule.params -gpu 1 -eggshell_triplet rays.txt

core.init: Mini-Rosetta version 49812:52965M from https://svn.rosettacommons.org/source/workspaces/johnk/rosetta/rosetta_source
core.init: command: /home/karenkhar/rosetta_source/bin/DARC.opencl.linuxgccrelease 
-database /home/karenkhar/rosetta_database/ -input_protein_file 2YXJ.pdb -input_ligand_file molecule.pdb 
-extra_res_fa molecule.params -gpu 1 -eggshell_triplet rays.txt
core.init: 'RNG device' seed mode, using '/dev/urandom', seed=-500381714 seed_offset=0 real_seed=-500381714
core.init.random: RandomGenerator:init: Normal mode, seed=-500381714 RG_type=mt19937
core.chemical.ResidueType: WARNING:: unrecognized residuetype property: BRANCH_POINT
core.chemical.ResidueType: WARNING:: unrecognized residuetype property: BRANCH_POINT
core.chemical.ResidueType: WARNING:: unrecognized residuetype property: BRANCH_POINT
core.chemical.ResidueType: WARNING:: unrecognized residuetype property: BRANCH_POINT
core.chemical.ResidueTypeSet: Finished initializing fa_standard residue type set.  Created 6413 residue types
core.conformation.Conformation: [ WARNING ] missing heavyatom:  OXT on residue ASN_p:CtermProteinFull 137
core.pack.task: Packer task: initialize from command line() 
basic.gpu: GPU #1: GeForce GTX 580 [1630 MHz] with 16 processors, 1024 threads.
core.pack.task: Packer task: initialize from command line() 
BEST FITNESS:	DARC_Zinc30_0_2YXJ_-1.pdb	14.2622
RMSD [Zinc30_0_2YXJ_-1]:26.1283

The score for the final docked pose is 14.2622. This is an indicator of how well the final docked pose sterically fits in the pocket.
The final docked pose is written to a pdb file name DARC_moleculeName_proteinName_-1.pdb
