-------------------------------------------------------------------------------------------------
   INSTRUCTIONS FOR CALCULATING EXPLICIT UNFOLDED STATE ENERGIES FOR NONCANONICAL AMINO ACIDS
-------------------------------------------------------------------------------------------------
KEYWORDS: NONCANONICALS ANALYSIS

Calculating the explicit unfolded state energies is the third of three steps toward being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles and HowToMakeRotamerLibraries. 

The explicit unfolded state energies of an amino acid represent the energy of an amino acid in the unfolded state of a protein and is used to replace the reference energies in Rosetta. The UnfoldedStateEnergyCalculator uses a fragment based method to calculate the average unfolded state energies for each ResidueType. The protocols works on a large set of protein structures that are split in to randomly generated fragments. The central residue of each fragment is mutated to the residue of interest. The fragment is repacked. The unweighted energy for each energy method in the scoring function is recorded for the central residue. After the energies for all fragment central residues are collected, a boltzmann-weighted-average average energy is calculated for each term. 

Calculation explicit unfolded state energies for a NCAA requires three steps:
 - Obtaining a set of input pdbs
 - Running the UnfoldedStateEnergyCalculator protocol on the set of pdbs
 - Modifying the unfolded state energies file in the database

-------------------------------------------
   STEP 01 OBTAINING A SET OF INPUT PDBS
-------------------------------------------

Since the UnfoldedStateEnergyCalculator protocol uses fragments from protein structures, we need a set of high quality structures to work with. Through their PISCES server, the Dunbrack laboratory maintains lists of structures in the Protein Data Bank organized based on xray resolution, precent sequence similarity, and r-factors**. These lists are a convenient way to get a set of high quality structures. Normally you would do the following:

1. Go the [PISCES server](http://dunbrack.fccc.edu/Guoli/pisces_download.php#culledpdb) and download e.g. the file that would show you pdbs with an xray resolution of at least 1.6 angstroms, less than 20% sequence identity, and r-factors of less than 0.25. You can use the provided script for that:

> To save space, there are a few pdbs provided. Feel free to delete them and download the hole list instead.

```
> cd inputs/
> rm *.pdb
> ../scripts/get_pdbs.bash cullpdb_pc20_res1.6_R0.25_d110520_chains1859
```
It will use the list to download many pdb files from the PDB. 

We use here also text file containing a list of them called cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned in the inputs directory. Rosetta will sometimes fail to correctly read in particular pdbs files. List the pdb files:

`ls inputs | grep '\.pdb' > list > inputs/list` 


**Citation: G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling server. Bioinformatics, 19:1589-1591, 2003. 

----------------------------------------------------------------
   STEP 02 RUNNING THE UNFOLDEDSTATEENERGYCALCULATOR PROTOCOL
----------------------------------------------------------------

The UnfoldedStateEnergyCalculator is relatively easy to run. The command line options are described bellow:

frag_size: single integer value, sets the number of residues in each fragment, should be an odd number and has a default of 5 which is what was used in the accompanying publication
residue_name: string value, sets the three letter code of the residue type which the central residue will be mutated to
repack_fragments: boolean value, controls if the fragments will be repacked before scoring and defaults to true
native_sequence: boolean value, controls if the central residue will be mutated before scoring and defaults to false

Additionally it is strongly recommended to add the following flags as they will make Rosetta handle more pdb files and improves runtime by disabling default features that will be negated by the fragmenting and prepacking

ignore_unrecognized_res: causes Rosetta to ignore unrecognized residue types and 
ex1 and ex2 and extrachi_cutoff 0: force rosetta to use additional rotamer during the fragment repacking
mute all and unmute devel.UnfoldedStateEnergyCalculator and unmute protocols.jd2.PDBJobInputer: reduces the size of the log file significantly by turning off unnecessary output
no_optH true: turns off the hydrogen optimization done when the protein is first read in 
detect_disulf false: turns off disulfide detection

Continuing the ornithine example we have used in the two previous protocol captures, to calculate the unfolded state energies one would run the following command.

```
$> $ROSETTA3/UnfoldedStateEnergyCalculator.linuxclangrelease -ignore_unrecognized_res -extrachi_cutoff 0 -l ../inputs/list -residue_name C40 -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_c40
```

NOTE: The extension on your executable my be different.

The run will take between 30-60 seconds per pdb file.

The log file contains lots of useful information. It contains the unweighted energies for each of the energy methods for each of the individual fragments. At the end it will print the average unweighted energies for each ResidueType as well as the Boltzmann weighted average unweighted energies. Boltzmann weighted average unweighted energies are used because some backbones just can't tolerate a mutation to a particular ResidueType and there are extremely high repulsive energies for some fragments that skew the average value. Using the Boltzmann weighting removes the higher energy outliers in a more elegant fashion than a hard energy cutoff.

-----------------------------------------------------
   STEP 03 MODIFY THE UNFOLDED STATE ENERGIES FILE
-----------------------------------------------------

Once the UnfoldedStateEnergyCalculator has finished running the Boltzmann weighted average unweighted energies need to be added to the database. The line you want is the "BOLZMANN UNFOLDED ENERGIES". These are the Boltzmann weighted average unfolded energies for each energy method. The file you need to modify is unfolded_state_residue_energies_mm_std.

Using the ornithine line as an example, the line form the log file is... 

BOLZMANN UNFOLDED ENERGIES:  fa_atr:    -2.462 fa_rep:     1.545 fa_sol:     1.166 mm_lj_intra_rep:     1.933 mm_lj_intra_atr:    -1.997 mm_twist:     2.733 pro_close:     0.009 hbond_sr_bb:    -0.006 hbond_lr_bb:     0.000 hbond_bb_sc:    -0.001 hbond_sc:     0.000 dslf_ss_dst:     0.000 dslf_cs_ang:     0.000 dslf_ss_dih:     0.000 dslf_ca_dih:     0.000

We could add the following to the unfolded_state_residue_energies_mm_std file in the database using the command bellow.

$ echo "C40 -2.462 1.545 1.166 1.933 -1.997 2.733 0.009 -0.006 0.000 -0.001  0.000" >> minirosetta_database/scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std 

The ResidueType can now be used in almost any Rosetta protocol that is compatible with the MM_STD scoring function.
---------------------------------------------------------------
##Part III -  INSTRUCTIONS FOR CREATING NONCANONICAL AMINO ACID ROTAMER LIBRARIES

KEYWORDS: NONCANONICALS UTILITIES

Creating a Noncanonical Amino Acid (NCAA) rotamer library is the second of two steps toward being able to use a NCAA in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles.  

Rotamer libraries are sets of common side chain conformations that generally correspond to local minima on the side chain conformational energy landscape. Side chain conformations are usually represented as a set of mean angles and a standard deviation to indicate variability. Rotamer libraries are used in for two main purposes in Rosetta: to provide starting points for side chain optimization routines, and the relative frequency is used as a pseudo-energy. Traditionally rotamer libraries are created by collecting statistics from protein structures. Rosetta uses the backbone dependent Drunbrack rotamer libraries. Since there are not enough structures containing NCAAs they must be generated.

Running the MakeRotLib protocol consists of four steps
 - creating and input template and generating the MakeRotLib options files
 - running the MakeRotLib protocol on each option file
 - Assembling the individual rotamer libraries in a single file
 - modify the ResidueType parameter file to be aware of our new rotamer library



###   STEP 01 MAKING INPUT FILES


Rosetta primarily uses backbone dependent rotamer libraries. Backbone-dependent rotamer libraries list provide side chain conformations sampled from residue positions whose backbone dihedral angles fall in particular bins. In the case of the Drunbrack rotamer libraries used by Rosetta the bins are in 10 degree intervals for for both phi and psi for a total of 1296 (36*36) phi/psi bins. To replicate this for the NCAAs we need to create a set of side chain rotamers for each member of a set of phi/psi bins.

The MakeRotLib protocol takes an option file as input. It requires an options file for each phi/psi bin. The first step in running it is creating these 1296 options files. Continuing from the HowToMakeResidueTypeParamFiles protocol capture we are again using ornithine as an example. Ornithine has 3 sidechain dihedral angles (chi). We want to sample each chi angle from 0 to 360 degrees in 30 degree intervals, and based on the chemistry of the side chain we predict that were will probably be three preferred angles for each chi angle at 60, 180, and 300 degrees for a total of 27 rotamers (3x3x3). We setup our MakeRotLib options file template as shown bellow.

```
<<<<< C40_rot_lib_options_XXX_YYY.in start >>>>>
AA_NAME C40
PHI_RANGE XXX XXX 0
PSI_RANGE YYY YYY 0
NUM_CHI 3
CHI_RANGE 1 0  330  30
CHI_RANGE 2 0  330  30
CHI_RANGE 3 0  330  30
CENTROID 300 1 300 1 300 1
CENTROID 300 1 300 1 180 2
CENTROID 300 1 300 1  60 3
CENTROID 300 1 180 2 300 1
CENTROID 300 1 180 2 180 2
CENTROID 300 1 180 2  60 3
CENTROID 300 1  60 3 300 1
CENTROID 300 1  60 3 180 2
CENTROID 300 1  60 3  60 3
CENTROID 180 2 300 1 300 1
CENTROID 180 2 300 1 180 2
CENTROID 180 2 300 1  60 3
CENTROID 180 2 180 2 300 1
CENTROID 180 2 180 2 180 2
CENTROID 180 2 180 2  60 3
CENTROID 180 2  60 3 300 1
CENTROID 180 2  60 3 180 2
CENTROID 180 2  60 3  60 3
CENTROID  60 3 300 1 300 1
CENTROID  60 3 300 1 180 2
CENTROID  60 3 300 1  60 3
CENTROID  60 3 180 2 300 1
CENTROID  60 3 180 2 180 2
CENTROID  60 3 180 2  60 3
CENTROID  60 3  60 3 300 1
CENTROID  60 3  60 3 180 2
CENTROID  60 3  60 3  60 3
<<<<< C40_rot_lib_options_XXX_YYY.in end >>>>>
```

AA_NAME <three letter code for the amno acid> 
PHI_RANGE <phi value for this bin> <phi value for this bin> 0 : The phi range functionality is not functional. Both values need to be the same and the interval set to 0
PSI_RANGE <psi value for this bin> <psi value for this bin> 0 : The psi range functionality is not functional. Both values need to be the same and the interval set to 0
NUM_CHI <number side chain dihedral angles> : This should be the same as in the parameter file.
CHI_RANGE <chi number> <starting value> <ending value> <interval> : The number of CHI_RANGE fields needs to equal the values specified for NUM_CHI.
CENTROID <Rotamer number for chi 1> <starting value> {<rotamer number for chi 2> <starting value>}{etc.} : CENTROIDS specify the starting points for the K-means clustering described in the related publication. A CENTROID field is needs for each potential rotamer. The number of CENTROID fields defines the number of rotamers listed in the resulting rotamer library.

To generate the 1296 input files we use a provided script that simply replace the XXX and YYY with the phi and psi values. The script is run as shown bellow.

```
$ cd inputs
$ ../scripts/make_inputs C40_rot_lib_options_XXX_YYY.in
```

The number of chi angles and the CHI_RANGE sampling interval are the primary determinants of the run time as they determine the number of rotamers that will be tested for each phi/psi bin. It is recommended to have at least 500 samples per chi. In the ornithine example we sample in 30 degree intervals for each of the 3 chi angles giving us a total of 1728 (12x12x12) conformations tested for each phi/psi bin. For a residue with a single chi 1 degree bins will suffice. 


###   STEP 02 RUNNING THE MAKEROTLIB PROTOCOL


The next step is to run the MakeRotLib protocol on each of the input files we created in step one. This is the most time consuming portion of the process and should probably be done on a cluster. As cluster setups vary, an example for a single MakeRotLib options file is provided. The other 1295 should be run identically.

```
$> ln -s HowToMakeRotamerLibraries/inputs .
$> $ROSETTA/bin/make_rot_lib.macosgccrelease  -rot_lib_options_file inputs/C40_rot_lib_options_-60_-40.in >& C40_rot_lib_options_-60_-40.log &
```

NOTE: The extension on your executable maybe different.

THe only options passed to the executable are the path to the database and the MakeRotLib options file. After the run completes a file called C40_-60_-40_bbdep.rotlib should be in the output directory. This is the backbone dependent rotamer library for a phi of -60 and a psi of -40.

The log file from the rosetta run in includes quite a bit of useful output.

There are three main sections to the log output,  "ROTAMERS", "CENTROIDS" and "FINAL ROTAMERS" sections. Each one shows the following data: phi, psi, omega and epsilon backbone dihedral angles, probability, total energy, torsion energy, intra-residue repulsive, intra-residue attractive, the number of chi angles, the assigned cluster number, the set of input chi angles, the set of minimized chi angles, the standard deviation, and the distance from that point to each of the cluster centroids.

The log file also displays the number of conformations per cluster, the average distance between the cluster center and the members of that cluster. Lack of conformations in a cluster and a large (>30) average cluster centroid distance suggests that that cluster is higher in energy. 


###   STEP 03 ASSEMBLING THE INDIVIDUAL ROTAMER LIBRARIES IN TO A SINGLE FILE


After the MakeRotLib protocol has been run on all of the MakeRotLib options files the individual rotamer libraries for each phi psi bin need to be assembled in to a single file. This is accomplished with a provided script as shown bellow. 

```
$> ln -s HowToMakeRotamerLibraries/scripts .
$> scripts/make_final_from_phi_psi.pl C40
```

The single file rotamer library should be called C40.rotlib. The file should be placed in the ncaa_rotlibs directory in the database. 

`$ cp C40.rotlibs Rosetta/main/database/ncaa_rotlibs/`


###   STEP 04 MODIFYING THE RESIDUE TYPE PARAM FILE


The last step is modifying the residue type parameter file to use the new rotamer library. To do this we need to add the name of the rotamer library, the number chi angles it describes, and how many bins there are for each chi angle to the Residue type parameter file. The ornithine rotamer library is called C40.rotlib and the rotamer library describes 3 chi angles and each of those 3 chi angle has 3 rotamer numbers. So we would use the following commands to add that information to the file we created in HowToMakeResidueTypeParamFiles.

```
$ echo "NCAA_ROTLIB_PATH C40.rotlib" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params

$ echo "NCAA_ROTLIB_NUM_ROTAMER_BINS 3 3 3 3" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params
```

