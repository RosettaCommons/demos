DARC Demo
=========

KEYWORDS: DOCKING LIGANDS 

Docking Approach using Ray Casting (DARC) is structure-based computational method for carrying out virtual screening by docking small-molecules into protein surface pockets.
In this demo, DARC is used to dock a small molecule in a pocket centered around residue 61 of the protein, E3 ubiquitin-protein ligase Mdm2 (PDB: 4ERF).

Generating input files
----------------------

* Protein files:

  4ERF.pdb was downloaded from The Protein Databank. Remove any water or ligand 
  or other HETATM present in the protein and remove redundant chains. Although 
  Rosetta will build any missing atoms including hydrogens in input protein PDB 
  files on the fly defore it uses them, we need to prepare the protein to 
  generate the electrostatic potential grid. This is one way to dump the 
  protein after building missing atom and added hydrogens. You can use the PDB provided in the input file to continue:
  
  **NOTE** (`$ROSETTA3`= path-to-Rosetta/main/source)
  
```
      $> cp input/4ERF.pdb .
      $> $ROSETTA3/bin/score.default.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false
```
  which gives the output protein 4ERF_0001.pdb

* Ligand files:

  Input files are included in the folder input, but the following explains how 
  they were obtained or generated. For more information, you can check the [prepare ligand tutorial](prepare_ligand) for more information. 
  
  4ERF.pdb was downloaded from The Protein 
  Databank. The small molecule used for docking was found in ZINC, the database 
  of commercially available compounds at http://zinc.docking.org/. We download 
  ZINC13989607 in mol2 format. The file name is `zinc_13989607.mol2`.

  Since the charges assigned to the atom is same for all conformers, it is fast 
  to assign charges before generating conformers. The charges were added to the 
  compound using the molcharge program from OpenEye (http://www.eyesopen.com/) 
  using the following command:

      $ OpenEye/bin/molcharge -in zinc_13989607.mol2 -out zinc_13989607_charged.mol2 -method am1bccsym

  If you have the compound in 2D SMILES format, then first convert to 3D format and then add charges.
  You can use openbabel or OpenEye program for converting from SMILES format to mol2 format

      $ OpenEye/bin/omega2 -in zinc_13989607_charged.mol2 -out zinc_13989607_conformers.mol2 -maxconfs 100

  then we need to generate parameter files for the compound to use in Rosetta.
  save the names of the compounds in a text file for generating parameters. You can use other softwares to perform this steps.
  
  Next, parameter files for the compounds are generated using batch_molfile_to_params.py which is a python app in rosetta. This command is used:

      $ echo zinc_13989607_conformers.mol2 > molfile_list.txt

      $ ROSETTA3/scripts/python/public/batch_molfile_to_params.py --script_path=<path-to-Rosetta>/main/source/scripts/python/public/molfile_to_params.py molfile_list.txt

  The expected output files are:

      params/zinc_13989607_conformers/000_conformers.pdb
      params/zinc_13989607_conformers/000.params

  You can copy them to your directory:
  ```
    $> cp input/000_conformers.pdb .
    $> cp input/000.params .
  ```

* Other input files:

  To run DARC we need to generate a RAY file for the input protein. To generate 
  this ray-file we need to input the protein in PDB format and specify a target 
  residue at the interface to generate the vectors. For more information about the ray file and how DARC works please refer to the provided citations. We can specify more than one residue at the 
  interface. The command to run DARC is based on shape only is as follows:
  **NOTE** For your special purposes, you will need to change the residue numbers and also provide params file for your ligands.

      $> $ROSETTA3/bin/make_ray_files.default.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only 

  Expected output from this command will be a ray-file named 
  `ray_4ERF_0001_61.txt`.

  To use a center of mass of any residue as origin points for casting rays:
```
      $> $ROSETTA3/bin/make_ray_files.default.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A
```
  Expected output: `ray_4ERF_0001_61.txt`

  To use multiple origin points for casting rays:
```
      $> $ROSETTA3/bin/make_ray_files.default.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61,54 -darc_shape_only -set_origin 5 -origin_res_num 85:A -multiple_origin
```
  Expected output: `ray_4ERF_0001_61,54.txt`

  To use a bound ligand to center the grid:
```
      $> $ROSETTA3/bin/make_ray_files.default.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -set_origin 5 -origin_res_num 85:A -bound_ligand input/4ERF_XTL_0001.pdb -extra_res_fa input/4ERF_XTL.params -lig_grid
```
  Expected output: `ray_4ERF_0001_61.txt`

  To include electrostatics calculations:

  For electrostatics calculations, we first need to generate the electrostatic potential grids. To generate the electrostatic potential grid, we can use the examples 
  Listings provided in the zap toolkit manual.

      $ OpenEye/bin/Listing_2 -in 4ERF_0001.pdb -out 4ERF.agd -buffer 2 -grid_spacing 0.5 -epsout 80 -epsin 1

  Expected output: `4ERF.agd`
  
  You can also copy this file from the input directory:
  ```
      $> cp input/4ERF.agd .
  ```
  
  Now we need to resize the electrostatic 
  potential grid (generated from openeye '4ERF.agd') to match the size of the 
  interface pocket grid. This step can be carried out while generating the ray 
  file.
```
      $> $ROSETTA3/bin/make_ray_files.linuxgccrelease -database Rosetta/main/database/ -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61,54 -set_origin 5 -origin_res_num 85:A -multiple_origin -bound_ligand input/4ERF_XTL_0001.pdb -extra_res_fa input/4ERF_XTL.params -lig_grid -espGrid_file 4ERF.agd
```
  The output from this command will be a ray-file named 
  `ray_4ERF_0001_61,54.txt` and an electrostatic potential grid file named 
  `DARC_4ERF.agd` which we use as input for running docking using DARC. 

Running DARC
------------

In the second step, we are actually running the docking calculations using the 
pre-generated ray-file. Call DARC for running the docking calculations using 
the pre-generated ray-file and the matched electrostatic potential grid. use 
the flag -num_particles to increase or decrease the sampling. high number of 
particles will take a longer time to finish. For best results, We suggest 
alleast 100 particles and 100 runs for docking single conformer or iterative 
docking (one-by-one) of ligand conformers. For on-the-fly sampling of ligand 
conformers we suggest at least 500 particles and 500 runs. 

Here we give theinput ligands for screening against the ray-file, as follows. Note that the options.short files are provided for faster test runs in case needed. 

Docking single ligand conformer:
```
    $>  $ROSETTA3/bin/DARC.default.linuxgccrelease @darc_single.options  
```
Expected output: `DARC_4ERF_0001_LG1.pdb`

Docking multiple ligand conformers:
```
    $>  $ROSETTA3/bin/DARC.linuxgccrelease @darc_conformers.options 
```
Expected output: `DARC_4ERF_0001_000.pdb`

To run DARC with shape only (without including electrostatics score):

    $>  $ROSETTA3/bin/DARC.linuxgccrelease @darc_shape.options

Expected output: `DARC_4ERF_0001_000.pdb`

To search conformers on-the-fly:

    $>  $ROSETTA3/bin/DARC.linuxgccrelease @darc_fly.options -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -search_conformers true

Expected output: `DARC_4ERF_0001_000.pdb`

The output of the DARC run is a docked model of the protein-ligand complex 
named `DARC_4ERF_0001_000.pdb`. To minimize the docked model: We can do the 
fullatom minimization of the DARC models separately in Rosetta or as an 
additional option while running DARC itself. To minimize the DARC models 
immediately after docking we add the flag `-minimize_output_complex`:

    $>  $ROSETTA3/bin/DARC.linuxgccrelease @darc_shape.options -minimize_output_complex

Expected output files are:

* DARC docked model complex: `DARC_4ERF_0001_000.pdb`
* DARC docked and minimized model complex: `mini_4ERF_0001_000.pdb`

Other options for DARC:

* `-use_connolly_surface`:  use this flag to make ray files with connolly surface insted of pocketshell generated by pocketgrid
* `-num_particles`:         set the number of particles used durin particle swarm optimization
* `-num_runs`:              set the number of particles used durin particle swarm optimization
* `-missing_point_weight`:  weight for rays misses the pocket
* `-steric_weight`:         weight for rays missesthe ligand
* `-extra_point_weight`:    weight for ligand moves away from the pocket
* `-esp_weight`:            electrostatics weight
* `-pocket_static_grid`:    no autoexpansion of the pocket grid (set ON for DARC)

Example: Commands for sample DARC run for protein MDM2 [PDB:4ERF]
-----------------------------------------------------------------

    $ cd run_dir

Copy 4ERF.pdb, 000_conformers.pdb, and 000.params, 4ERF.agd, DARC_4ERF.agd, 4ERF_0001.pdb, 4ERF_XTL_0001.pdb, 4ERF_XTL.params to run_dir. 
Run DARC with the following command:

* DARC shape only:

        $ Rosetta/main/source/bin/score.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false
        $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only
        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -darc_shape_only -num_particles 50 -num_runs 50 

  Expected output: `DARC_4ERF_0001_LG1.pdb`

* DARC shape and electrostatics:

        $ Rosetta/main/source/bin/make_ray_files.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -espGrid_file 4ERF.agd
        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 50 -num_runs 50

  Expected output: `DARC_4ERF_0001_LG1.pdb`

* DARC sampling conformers on the fly

        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 500 -num_runs 500 -search conformers true

  Expected output: `DARC_4ERF_0001_000.pdb`

* DARC sampling conformers iteratively one-by-one

        $ Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 500 -num_runs 500 -search conformers false

  Expected output: `DARC_4ERF_0001_000.pdb`

Running DARC in GPU systems
---------------------------
DARC is adapted to run on GPU systems. On average, Running DARC on GPU system 
is 190 fold faster than running on CPU. The results from DARC that runs on CPU 
and GPU systems should be same (with same input and constant seed).

First make sure you built rosetta with the flag `extras=opencl` to run on GPU 
systems

    $ cd Rosetta/main/source
    $ ./scons.py mode=release extras=opencl bin


Then use the opencl executables (eg: 
DARC.opencl.linuxgccrelease) instead of default executables with the flag '-gpu 1', i.e call DARC.opencl.linuxgccrelease 
instead of DARC.linuxgccrelease. All other options are same as CPU.
Here is an example for GPU command:

    $ Rosetta/main/source/bin/score.opencl.linuxgccrelease -in:file:s 4ERF.pdb -out:output -no_optH false -gpu 1
    $ Rosetta/main/source/bin/make_ray_files.opencl.linuxgccrelease -pocket_static_grid -protein 4ERF_0001.pdb -central_relax_pdb_num 61 -darc_shape_only -gpu 1
    $ Rosetta/main/source/bin/DARC.opencl.linuxgccrelease -protein 4ERF_0001.pdb -ligand 000_conformers.pdb -extra_res_fa 000.params -ray_file ray_4ERF_0001_61.txt -darc_shape_only -gpu 1 -num_particles 500 -num_runs 500 -search_conformers true

Expected output: `DARC_4ERF_0001_000.pdb`

Running DARC will result in these output files:

* DARC_4ERF_0001_LG1.pdb is the final docked pose for single conformer.
* DARC_4ERF_0001_000.pdb is the final docked pose for multiple conformers. 
  Note:The flag '-minimize_output_complex' can be used to minimized the output 
  model
* mini_4ERF_0001_LG1.pdb is final docked and minimized pose for single 
  conformers.
* mini_4ERF_0001_000.pdb is final docked and minimized pose for multiple 
  conformers.

Analyzing the results
---------------------

The output should look similar to:

    Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 10 -num_runs 10
    core.init: Rosetta version 1e0e3a1462814b2fa93580e9e1b1f7ab81736efc 2014-11-19 15:22:07 -0600 from git@github.com:RosettaCommons/main.git
    core.init: command: /Users/ragul/Rosetta/main/source/bin/DARC.linuxgccrelease -protein 4ERF_0001.pdb -ligand 4ERF_XTL_0001.pdb -extra_res_fa 4ERF_XTL.params -ray_file ray_4ERF_0001_61.txt -espGrid_file DARC_4ERF.agd -num_particles 10 -num_runs 10
    core.init: 'RNG device' seed mode, using '/dev/urandom', seed=590325598 seed_offset=0 real_seed=590325598
    core.init.random: RandomGenerator:init: Normal mode, seed=590325598 RG_type=mt19937
    core.init: Resolved executable path: /Users/ragul/Rosetta/main/source/bin/DARC.linuxgccrelease
    core.init: Looking for database based on location of executable: /Users/ragul/Rosetta/main/database/
    core.chemical.ResidueTypeSet: Finished initializing fa_standard residue type set.  Created 738 residue types
    core.pack.task: Packer task: initialize from command line() 
    origin_space : 7.25
    Reading ligand 4ERF_XTL_0001.pdb
    core.pack.task: Packer task: initialize from command line() 
    Starting PSO
    Time for DARC runs: 0.589119
    SCORES : 4ERF_0001_LG1	Conformer 1 has DARC Score : 4.29993

In this case the DARC score for the final docked pose is 4.29993. This score 
can be use to compare with other ligand that are docked using the same ray 
file. The lower the score the better the match between the protein pocket and 
ligand. `DARC_4ERF_0001_LG1.pdb` is the out file that contains the final 
docked pose for single conformer.

If we use the flag `-minimize_output_complex` the model will be minimized and 
we get the file `mini_4ERF_0001_LG1.pdb` as output.

