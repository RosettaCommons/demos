--------------------------------------------------------------------------------
   INSTRUCTIONS FOR RUNNING DARC PROTOCOL
--------------------------------------------------------------------------------

	Docking Approach using Ray Casting (DARC) is structure-based computational method 
for carrying out virtual screening by docking small-molecules into protein surface pockets.
IN this demo, DARC is used to dock a small molecule in a pocket centered around residue 61 
of the protein, E3 ubiquitin-protein ligase Mdm2  (PDB:4ERF). 

--------------------------------------------------------------------------------
   GENERATING INPUT FILES
--------------------------------------------------------------------------------

PROTEIN FILES:

4ERF.pdb was downloaded from The Protein Databank. Remove any water or lignad or other HETATM 
present in the protein and remove redundant chains. Although Rosetta will build any missing atoms 
including hydrogens in input protein PDB files on the fly defore it uses them, we need to prepare 
the protein to generate the electrostatic potential grid. 
This is one way to dump the protein after building missing atoma and added hydrogens:
$ Rosetta/main/source/bin/score.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false
which gives the output protein 4ERF_0001.pdb

LIGAND FILES:

	Input files are included in the folder inputs, but the following explains 
how they were obtained or generated. 4ERF.pdb was downloaded from The Protein Databank.
The small molecule used for docking was found in ZINC, the database of commercially 
available compounds at http://zinc.docking.org/. We download ZINC13989607 in mol2 format. 
The file name is zinc_13989607.mol2

Since the charges assigned to the atom is same for all conformers, it is fast to assign charges before generating conformers. 
The charges were added to the compound using the molcharge program from OpenEye 
http://www.eyesopen.com/ using the following command:

$ OpenEye/bin/molcharge -in zinc_13989607.mol2 -out zinc_13989607_charged.mol2 -method am1bccsym

If you have the compound in 2D SMILES format, then first convert to 3D format and then add charges.
You can use openbabel or OpenEye program for converting from SMILES format to mol2 format

$ OpenEye/bin/omega2 -in zinc_13989607_charged.mol2 -out zinc_13989607_conformers.mol2 -maxconfs 100

then we need to generate parameter files for the compound to use in Rosetta.
save the names of the compounds in a text file for generating parameters. Next, parameter files for the compounds are generated using batch_molfile_to_params.py which is a python app in rosetta. This command is used:For Ex:

$ echo zinc_13989607_conformers.mol2 > molfile_list.txt

$ Rosetta/main/source/src/python/apps/public/batch_molfile_to_params.py -d Rosetta/main/database --script_path=Rosetta/main/source/src/python/apps/public/molfile_to_params.py molfile_list.txt

The output files are
params/zinc_13989607_conformers/000_conformers.pdb
params/zinc_13989607_conformers/000.params

OTHER INPUT files for DARC:

To run DARC we need to generate a RAY file for the input protein. To generate this ray-file we need to input the protein in PDB format and specify a target residue at the interface. We can specify more than one residue at the interface. The command to run DARC is as follows:

$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only 

The output from this command will be a ray-file named “ray_4ERF_0001_61.txt”.

To use a center of mass of any residue as origin points for casting rays:
$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A
output : ray_4ERF_0001_61.txt

To use multiple origin points for casting rays:

$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A -multiple_origin
output : ray_4ERF_0001_61,54.txt

To use a bound ligand to center the grid:

$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A -bound_ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -lig_grid
output : ray_4ERF_0001_61.txt

To include electrostatics calculations:

For electrostatics calculations, first we need to resize the electrostatic potential grid (generated from openeye '4ERF.agd') to match the size of the interface pocket grid. This step can be carried out while generating the ray file.

To generate the electrostatic potential grid, we can use the examples Listings provided in the zap toolkit manual.
$ OpenEye/bin/Listing_2 -in 4ERF_0001.pdb -out 4ERF.agd -buffer 2 -grid_spacing 0.5 -epsout 80 -epsin 1
output : 4ERF.agd

$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61,54 -set_origin 5 -origin_res_num 85:A -multiple_origin -bound_ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -lig_grid -espGrid_file 4ERF.agd
output : ray_4ERF_0001_61,54.txt
output : DARC_4ERF.agd

The output from this command will be a ray-file named “ray_4ERF_0001_61,54.txt” and an electrostatic potential grid file named “DARC_4ERF.agd” which we use as input for running docking using DARC. 
 
--------------------------------------------------------------------------------
   RUNNING DARC
--------------------------------------------------------------------------------

In the second step, we are actually running the docking calculations using the pre-generated ray-file. Call DARC for running the docking calculations using the pre-generated ray-file and the matched electrostatic potential grid. use the flag -num_particles to increase or decrease the sampling. high number of particles will take a longer time to finish. 
For best results, We suggest alleast 100 particles and 100 runs for docking single conformer or iterative docking (one-by-one) of ligand conformers. For on-the-fly sampling of ligand conformers we suggest atleast 500 particles and 500 runs. Here we give the input ligands for screening against the ray-file, as follows: 

dockcing single ligand conformer:
$  Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF_0001.agd 

dockcing multiple ligand conformers:
$  Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF_0001.agd 

TO run DARC with shape only (without including electrostatics score):

$  Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file -darc_shape_only

TO search conformers on-the-fly:
$  Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF_0001.agd -search_conformers true

The output of the DARC run is a docked model of the protein-ligand complex named “DARC_4ERF_0001_000.pdb”

To minimize the docked model:
We can do the fullatom minimization of the DARC models separately in Rosetta or as an additional option while running DARC itself. To minimize the DARC models immediately after docking we add the flag “-minimize_output_complex”:

$  Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF_0001.agd -minimize_output_complex

--------------------------------------------------------------------------------
   OTHER OPTIONS FOR DARC
--------------------------------------------------------------------------------
-use_connolly_surface   : use this flag to make ray files with connolly surface insted of pocketshell generated by pocketgrid
-num_particles          : set the number of particles used durin particle swarm optimization
-num_runs		: set the number of particles used durin particle swarm optimization
-missing_point_weight	: weight for rays misses the pocket
-steric_weight		: weight for rays missesthe ligand
-extra_point_weight	: weight for ligand moves away from the pocket
-esp_weight 		: electrostatics weight
-pocket_static_grid 	: no autoexpansion of the pocket grid (set ON for DARC)
 
--------------------------------------------------------------------------------
   BUILD ROSETTA
--------------------------------------------------------------------------------

$ scons mode=release bin

To build with GPU enabled:
$ scons mode=release extras=opencl bin

--------------------------------------------------------------------------------
   EXAMPLE: COMMANDS FOR SAMPLE DARC RUN FOR PROTEIN MDM2 [PDB:4ERF]
--------------------------------------------------------------------------------

$ cd run_dir	

Copy 4ERF.pdb, 000_conformers.pdb, and 000.params, 4ERF.agd, DARC_4ERF.agd, 4ERF_0001.pdb, 4ERF_XTL_0001.pdb, 4ERF_XTL.params to run_dir. 
Run DARC with the following command:

DARC SHAPE ONLY:

$ Rosetta/main/source/bin/score.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false
$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only
$ Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -darc_shape_only -num_particles 50 -num_runs 50 

DARC SHAPE AND ELECTROSTATICS:

$ Rosetta/main/source/bin/make_ray_files.macosgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -espGrid_file 4ERF.agd
$ Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 50 -num_runs 50

DARC SAMPLING CONFORMERS ON THE FLY

$ Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 500 -num_runs 500 -search conformers true

DARC SAMPLING CONFORMERS ITERATIVELY ONE-BY-ONE

$ Rosetta/main/source/bin/DARC.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 500 -num_runs 500 -search conformers false

--------------------------------------------------------------------------------
  RUN DARC IN GPU SYSTEMS
--------------------------------------------------------------------------------

DARC is adapted to run on GPU systems. On average, Running DARC on GPU system is 190 fold faster than running on CPU.
The reuslts from DARC that runs on CPU and GPU systems should be same (with same input and constant seed). 

First make sure you built rosetta with the flag 'extras=opencl' to run on GPU systems

$ cd Rosetta/main/source
$ ./scons.py mode=release extras=opencl bin

Then use the opencl executables instead of default executables with the flag '-gpu 1', i.e call DARC.opencl.macosgccrelease 
instead of DARC.macosgccrelease. All other options are same as CPU.
Here is an example for GPU command:

$ Rosetta/main/source/bin/score.opencl.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false -gpu 1
$ Rosetta/main/source/bin/make_ray_files.opencl.macosgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -gpu 1
$ Rosetta/main/source/bin/DARC.opencl.macosgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -darc_shape_only -gpu 1 -num_particles 500 -num_runs 500 -search_conformers true


Output files:

Running DARC will result in this output file:
DARC_4ERF_0001_LG1.pdb is the final docked pose for single conformer.
DARC_4ERF_0001_000.pdb is the final docked pose for multiple conformers.


--------------------------------------------------------------------------------
   ANALYZING THE RESULTS
--------------------------------------------------------------------------------

The output should look similar to:

Rosetta/main/source/bin/DARC.macosclangrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 10 -num_runs 10
core.init: Rosetta version 1e0e3a1462814b2fa93580e9e1b1f7ab81736efc 2014-11-19 15:22:07 -0600 from git@github.com:RosettaCommons/main.git
core.init: command: /Users/ragul/Rosetta/main/source/bin/DARC.macosclangrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 10 -num_runs 10
core.init: 'RNG device' seed mode, using '/dev/urandom', seed=590325598 seed_offset=0 real_seed=590325598
core.init.random: RandomGenerator:init: Normal mode, seed=590325598 RG_type=mt19937
core.init: Resolved executable path: /Users/ragul/Rosetta/main/source/bin/DARC.macosclangrelease
core.init: Looking for database based on location of executable: /Users/ragul/Rosetta/main/database/
core.chemical.ResidueTypeSet: Finished initializing fa_standard residue type set.  Created 738 residue types
core.pack.task: Packer task: initialize from command line() 
origin_space : 7.25
Reading ligand 4ERF_XTL_0001.pdb
core.pack.task: Packer task: initialize from command line() 
Starting PSO
Time for DARC runs: 0.589119
SCORES : 4ERF_0001_LG1	Conformer 1 has DARC Score : 4.29993

In this case the DARC score for the final docked pose is 4.29993. This score can be use to compare with other ligand that are docked using the same ray file. The lower the score the better the match between the protein pocket and ligand. 'DARC_4ERF_0001_LG1.pdb' is the out file that containfs the final docked pose for single conformer.

If we use the flag '-minimize_output_complex' the model will be minimized and we get the file 'mini_4ERF_0001_LG1.pdb' as output.

