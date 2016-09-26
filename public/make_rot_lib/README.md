Making non-canonical amino acid rotamer libraries
=================================================

KEYWORDS: NONCANONICALS UTILITIES

This demo was written by P. Douglas Renfrew (renfrew at nyu dot edu).

This demo illustrates the make_rot_lib protocol. It was originally written as 
the 2nd (of 4) parts of a "protocol capture" acompanying the paper "Using 
Noncanonical Amino Acids in Computational Protein-Peptide Interface Design" by 
P. Douglas Renfrew, Eun Jung Choi, Brian Kuhlman (2011) PLoS One, from the 
RosettaCon2010 special PLoS ONE special issue. The four separate protocol 
captures that describe the four different ways Rosetta was used in the 
publication: creating Rosetta ResidueType parameter files for NCAAs, creating 
backbone dependent rotamer libraries for NCAAs, calculating explicit unfolded 
state reference energies for NCAAs, and the running of the 
DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol caputres 
for creating ResidueType paramater files, rotamer libraries and explicite 
unfolded state energies describe the process used in the publication but are 
written from the standpoint of a researcher looking to add an additional NCAA 
to Rosetta. 

Creating a Noncanonical Amino Acid (NCAA) rotamer library is the second of two 
steps toward being able to use a NCAA in Rosetta. To add a new NCAA or to 
better understand how the NCAAs in the related publication were added one 
should have already completed or understand the steps in 
HowToMakeResidueTypeParamFiles.

Rotamer libraries are sets of common side chain conformations that generally 
correspond to local minima on the side chain conformational energy landscape. 
Side chain conformations are usually represented as a set of mean angles and a 
standard deviation to indicate variability. Rotamer libraries are used in for 
two main purposes in Rosetta: to provide starting points for side chain 
optimization routines, and the relative frequency is used as a pseudo-energy. 
Traditionally rotamer libraries are created by collecting statistics from 
protein structures. Rosetta uses the backbone dependent Drunbrack rotamer 
libraries. Since there are not enough structures containing NCAAs they must be 
generated.

Running the MakeRotLib protocol consists of four steps:

1. Creating and input template and generating the MakeRotLib options files.
2. Running the MakeRotLib protocol on each option file.
3. Assembling the individual rotamer libraries in a single file.
4. Modifying the ResidueType parameter file to be aware of our new rotamer 
   library.

Step 1: Making input files
--------------------------

Rosetta primarily uses backbone dependent rotamer libraries. Backbone-dependent 
rotamer libraries list provide side chain conformations sampled from residue 
positions whose backbone dihedral angles fall in particular bins. In the case 
of the Drunbrack rotamer libraries used by Rosetta the bins are in 10 degree 
intervals for for both phi and psi for a total of 1296 (36x36) phi/psi bins. 
To replicate this for the NCAAs we need to create a set of side chain rotamers 
for each member of a set of phi/psi bins.

The MakeRotLib protocol takes an option file as input. It requires an options 
file for each phi/psi bin. The first step in running it is creating these 1296 
options files. Continuing from the HowToMakeResidueTypeParamFiles protocol 
capture we are again using ornithine as an example. Ornithine has 3 sidechain 
dihedral angles (chi). We want to sample each chi angle from 0 to 360 degrees 
in 30 degree intervals, and based on the chemistry of the side chain we predict 
that were will probably be three preferred angles for each chi angle at 60, 
180, and 300 degrees for a total of 27 rotamers (3x3x3). We setup our 
MakeRotLib options file template as shown bellow.

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

Field definitions:

* AA_NAME \<three letter code for the amno acid> 

* PHI_RANGE \<phi value for this bin> \<phi value for this bin> 0 : The phi range 
  functionality is not functional. Both values need to be the same and the 
  interval set to 0

* PSI_RANGE \<psi value for this bin> \<psi value for this bin> 0 : The psi range 
  functionality is not functional. Both values need to be the same and the 
  interval set to 0

* NUM_CHI \<number side chain dihedral angles> : This should be the same as in 
  the parameter file.

* CHI_RANGE \<chi number> \<starting value> \<ending value> \<interval> : The 
  number of CHI_RANGE fields needs to equal the values specified for NUM_CHI.

* CENTROID \<Rotamer number for chi 1> \<starting value> {\<rotamer number for chi 
  2> \<starting value>}{etc.} : CENTROIDS specify the starting points for the 
  K-means clustering described in the related publication. A CENTROID field is 
  needs for each potential rotamer. The number of CENTROID fields defines the 
  number of rotamers listed in the resulting rotamer library.

To generate the 1296 input files we use a provided script that simply replace 
the XXX and YYY with the phi and psi values. The script is run as shown bellow.

    $> cd inputs
    $> ../scripts/make_inputs.bash C40_rot_lib_options_XXX_YYY.in

The number of chi angles and the CHI_RANGE sampling interval are the primary 
determinants of the run time as they determine the number of rotamers that will 
be tested for each phi/psi bin. It is recommended to have at least 500 samples 
per chi. In the ornithine example we sample in 30 degree intervals for each of 
the 3 chi angles giving us a total of 1728 (12x12x12) conformations tested for 
each phi/psi bin. For a residue with a single chi 1 degree bins will suffice. 

Step 2: Running the MakeRotLib protocol
----------------------------------------

The next step is to run the MakeRotLib protocol on each of the input files we 
created in step one. This is the most time consuming portion of the process and 
should probably be done on a cluster. As cluster setups vary, an example for a 
single MakeRotLib options file is provided. The other 1295 should be run 
identically.

    $> cd outputs
    $> PATH/TO/ROSETTA/bin/MakeRotLib.macosgccrelease -rot_lib_options_file ../inputs/C40_rot_lib_options_-60_-40.in >& C40_rot_lib_options_-60_-40.log &

NOTE: The extension on your executable maybe different.

The only options passed to the executable are the path to the database and the 
MakeRotLib options file. After the run completes a file called 
C40_-60_-40_bbdep.rotlib should be in the output directory. This is the 
backbone dependent rotamer library for a phi of -60 and a psi of -40.

The log file from the rosetta run in includes quite a bit of useful output. 
There are three main sections to the log output,  "ROTAMERS", "CENTROIDS" and 
"FINAL ROTAMERS" sections. Each one shows the following data: phi, psi, omega 
and epsilon backbone dihedral angles, probability, total energy, torsion 
energy, intra-residue repulsive, intra-residue attractive, the number of chi 
angles, the assigned cluster number, the set of input chi angles, the set of 
minimized chi angles, the standard deviation, and the distance from that point 
to each of the cluster centroids. The log file also displays the number of 
conformations per cluster, the average distance between the cluster center and 
the members of that cluster. Lack of conformations in a cluster and a large 
(>30) average cluster centroid distance suggests that that cluster is higher in 
energy. 

Step 3: Assembling the individual rotamer libraries in to a single file
-----------------------------------------------------------------------

After the MakeRotLib protocol has been run on all of the MakeRotLib options 
files the individual rotamer libraries for each phi psi bin need to be 
assembled in to a single file. This is accomplished with a provided script as 
shown bellow. 

    $> cd outputs
    $> ../scripts/make_final_from_phi_psi.pl C40

The single file rotamer library should be called C40.rotlib. The file should be 
placed in the ncaa_rotlibs directory in the database. 

    > cp C40.rotlibs PATH/TO/DATABASE/rosetta_database/ncaa_rotlibs/

Step 4: Modifying the residue type PARAMS file
-----------------------------------------------

The last step is modifying the residue type parameter file to use the new 
rotamer library. To do this we need to add the name of the rotamer library, the 
number chi angles it describes, and how many bins there are for each chi angle 
to the Residue type parameter file. The ornithine rotamer library is called 
C40.rotlib and the rotamer library describes 3 chi angles and each of those 3 
chi angle has 3 rotamer numbers. So we would use the following commands to add 
that information to the file we created in HowToMakeResidueTypeParamFiles.

    > echo "NCAA_ROTLIB_PATH C40.rotlib" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params
    > echo "NCAA_ROTLIB_NUM_ROTAMER_BINS 3 3 3 3" >> PATH/TO/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/ornithine.params

