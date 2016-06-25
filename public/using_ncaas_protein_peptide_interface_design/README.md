# Using NCAAs in Protein-Peptide Interface Design

KEYWORDS: NONCANONICALS DOCKING PEPTIDES INTERFACES

## Author(s)
- P. Douglas Renfrew (renfrew@nyu.edu)
- Eun Jung Choi
- Brian Kuhlman

## Reference
[[Incorporation of Noncanonical Amino Acids into Rosetta and Use in Computational Protein-Peptide Interface Design|http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032637]]

## Brief Description
This is four separate protocol captures that describe the four different ways Rosetta was used in the accompanying publication: creating Rosetta ResidueType parameter files for NCAAs, creating backbone dependent rotamer libraries for NCAAs, calculating explicit unfolded state reference energies for NCAAs, and the running of the DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol caputres for creating ResidueType paramater files, rotamer libraries and explicite unfolded state energies describe the process used in the publication but are written from the standpoint of a researcher looking to add an additional NCAA to Rosetta. 

## Running
See README files for individual files for descriptions about each part was run. 

## Version
The code for this protocol capture is completely checked into trunk rXXXXX and database rYYYYY. 

-



## Part I -  CREATING ROSETTA RESIDUETYPE PARAMETER FILES FOR NONCANONICAL AMINO ACIDS

KEYWORDS: NONCANONICALS UTILITIES
Creating a Residue Type Parameter file is the first of three steps to being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand the how the NCAAs in the accompanying publication were added one should follow the steps in this folder followed by the steps in HowToMakeNCAARotamerLibraries.  

ResidueType Parameter files are the basic input file for the creation the core::chemical::ResidueType in Rosetta. ResidueTypes are used for storing information about the atoms, bonds, and chemical connectivity of both polymer (DNA/RNA, peptide) and ligand residues in Rosetta. Including an idealized conformation in an internal coordinate format. They also contain parameters used by some of the energy methods, and rotamer libraries.

This directory includes modified open babble source, molfile2params_polymer.py script, instructions, inputs and outputs for creating a param file for a alpha-peptide backbone with a noncanonical side chain. The example alpha-peptide is ornithine.

There are five basics steps to generating a ResidueType parameter file as in the accompanying publication outlined bellow:
- Creating an initial structure
- Minimizing that structure
- Converting it to a hybrid molfile 2000/3000 format
- Modifying the molfile with additional information used by the molfile_2_params_polymer.py script
- Running the molfile_2_params_polymer.py script

-


###   STEP 00 BUILDING OPENBABLE

Openbabble is an open source program distributed under the GLPv2 license used to create a specially formated molfile to be used as input for the molfile_2_params_plymer.py script. The copy included here has been modified to produce molfile 2000 format files with the molfile 3000 bond types. 

You can build the modified open babble source using the standard configure, make, make install. I like to install it locally but you can install it anywhere.

Note: The path for the configure script will need to change on your system.

```
$ cd openbable
$ mkdir install
$ tar xvzf openbabel-2.2.0.tar.gz
$ cd openbabel-2.2.0
$ ./configure --prefix="/PATH/TO/HowToMakeResidueTypeParamFiles/openbabel/install/" 
$ make
$ make install
```
the executable babble should now be in /PATH/TO/HowToMakeResidueTypeParamFiles/openbabel/install/bin/

------------------------------------------
###   STEP 01 GENERATING INITIAL STRUCTURE

The first step is to generate an initial structure. I usually do this in PyMOL but any molecular editor will work. PyMOL has rudimentary but very functional editing capabilities. The goal is to get the structure close to the final structure. Don't stress too much because the minimization will clean up most things. Make sure you double check chirality. To make sure the geometry for the atoms near the ends of the molecule is correct you will need to put capping groups on the molecule. I would recommend an acetyl (ACE) on the nitro-terminus of and a n-methyl (NME) on the carboxy terminus. This structure (ACE-X-NME) is usually called a dipeptide despite the fact that it only has one complete residue.

PRO TIPS:
 - Within PyMOL go to the "Mouse" drop down and select "3 button editing".
 - Within PyMOL type "help edit_keys" in the command box will bring up the instructions for the editor. 
 - Within PyMOL type "set valence, 1" in the command box to show single, double, triple bonds.
 - If you are making multiple ResidueTypes, be consistent with atom order. It makes future steps easier. I like to put all the atoms for the capping groups first, followed by backbone heavy atoms, side chain heavy atoms, and then then hydrogens.

EXAMPLE:
See the example pdb files in the folder stage_01_initial_structures and shown bellow...

```
<<<<< ornithine.pdb start >>>>>
ATOM      1  C   ACE     1       0.984   0.045  -0.578  1.00  0.00           C
ATOM      2  O   ACE     1       1.815  -0.721  -1.083  1.00  0.00           O
ATOM      3  CH3 ACE     1      -0.445  -0.396  -0.349  1.00  0.00           C
ATOM      4 1HH3 ACE     1      -0.468  -1.251   0.313  1.00  0.00           H
ATOM      5 2HH3 ACE     1      -1.013   0.407   0.098  1.00  0.00           H
ATOM      6 3HH3 ACE     1      -0.904  -0.668  -1.288  1.00  0.00           H
ATOM     26  N   NME     3       4.937   2.266   0.684  1.00  0.00           N
ATOM     27  CH3 NME     3       5.341   3.324   1.592  1.00  0.00           C
ATOM     28  H   NME     3       5.534   1.511   0.316  1.00  0.00           H
ATOM     29 1HH3 NME     3       4.689   4.177   1.478  1.00  0.00           H
ATOM     30 2HH3 NME     3       5.286   2.977   2.614  1.00  0.00           H
ATOM     31 3HH3 NME     3       6.355   3.627   1.379  1.00  0.00           H
ATOM      7  N   ALA     2       1.392   1.405  -0.195  1.00  0.00           N
ATOM      8  CA  ALA     2       2.806   1.518  -0.539  1.00  0.00           C
ATOM      9  C   ALA     2       3.512   2.478   0.390  1.00  0.00           C
ATOM     10  O   ALA     2       2.929   3.435   0.912  1.00  0.00           O
ATOM     11  CB  ALA     2       2.895   1.945  -2.014  1.00  0.00           C
ATOM     12  C01 ALA     2       2.214   0.926  -2.946  1.00  0.00           C
ATOM     13  N01 ALA     2       1.680   0.434  -5.295  1.00  0.00           N
ATOM     14  C02 ALA     2       2.334   1.404  -4.405  1.00  0.00           C
ATOM     15  H   ALA     2       0.770   2.167   0.251  1.00  0.00           H
ATOM     16  HA  ALA     2       3.283   0.529  -0.413  1.00  0.00           H
ATOM     17 2HB  ALA     2       2.409   2.923  -2.195  1.00  0.00           H
ATOM     18 3HB  ALA     2       3.943   2.035  -2.355  1.00  0.00           H
ATOM     19  H01 ALA     2       2.698  -0.045  -2.840  1.00  0.00           H
ATOM     20  H02 ALA     2       1.461  -0.405  -4.777  1.00  0.00           H
ATOM     21  H03 ALA     2       1.161   0.837  -2.679  1.00  0.00           H
ATOM     22  H04 ALA     2       2.300   0.206  -6.059  1.00  0.00           H
ATOM     23  H05 ALA     2       3.387   1.492  -4.674  1.00  0.00           H
ATOM     24  H06 ALA     2       1.850   2.375  -4.509  1.00  0.00           H
ATOM     25  H07 ALA     2       0.828   0.835  -5.661  1.00  0.00           H
END
<<<<< ornithine.pdb end >>>>>
```
-----------------------------------
###   STEP 02 MAKING GAUSSIAN INPUT

Rosetta assumes that We need to minimize the initial structure we made in PyMOL to get a good set of ideal bond lengths and angles. We will use Gaussian, an is a ab initio quantum mechanics package from Schodinger Inc., to do this but You can use the molecular modeling program of your choice. To do this we will take the coordinates from the pdb file and put them into a gaussian input file. Like the one shown bellow and in stage_02_gaussian_input. A discussion of the complete gaussian input structure is beyond the scope of this document. Documentation for gaussian can be found here http://www.gaussian.com/ . 

In short the input is as follows...
line 1 sets the path for the checkpoint file
line 2 describes the level of theory, options, and convergence criteria. You may need to change the basis set to something smaller if you are using stuff bellow the 4th line of the periodic table. 
lines 3-5 are comments
line 6 is the charge and multiplicity (this is usually 0 1 but ornithine is charged)
line 7-37 are the elemental type and xyz coordinates of the atoms from the pdb file in the same order as the pdb file
line 38 is blank
line 39-40 is the modredundant input that says that we want to keep the torsion formed by atoms 1 13 14 and 15 fixed at 150.00 degrees
line 41 is blank.

Running Gaussian is simple but the minimizations take a long time (a few hours per structure). The command bellow will run gaussian on the input file in the stage_02 folder and put the output in the stage_03 folder.

```
> g03 stage_02_gaussian_input/ornithine.com stage_03_gaussian_output/ornithine.log
```

```
<<<<< ornithine.com start >>>>>
%Chk=stage02_ace_nme_res_ordered_pdbs/ornithine.chk
\# HF/6-31G(d) Opt=ModRedundant SCF=Tight Test 

scan rotamers
 
1  1
C        0.984   0.045  -0.578
O        1.815  -0.721  -1.083
C       -0.445  -0.396  -0.349
H       -0.468  -1.251   0.313
H       -1.013   0.407   0.098
H       -0.904  -0.668  -1.288
N        4.937   2.266   0.684
C        5.341   3.324   1.592
H        5.534   1.511   0.316
H        4.689   4.177   1.478
H        5.286   2.977   2.614
H        6.355   3.627   1.379
N        1.392   1.405  -0.195
C        2.806   1.518  -0.539
C        3.512   2.478   0.390
O        2.929   3.435   0.912
C        2.895   1.945  -2.014
C        2.214   0.926  -2.946
N        1.680   0.434  -5.295
C        2.334   1.404  -4.405
H        0.770   2.167   0.251
H        3.283   0.529  -0.413
H        2.409   2.923  -2.195
H        3.943   2.035  -2.355
H        2.698  -0.045  -2.840
H        1.461  -0.405  -4.777
H        1.161   0.837  -2.679
H        2.300   0.206  -6.059
H        3.387   1.492  -4.674
H        1.850   2.375  -4.509
H        0.828   0.835  -5.661

1 13 14 15  -150.00 F
13 14 15 7  150.00 F

<<<<< ornithine.com end >>>>>
```
---------------------------------
###   STEP 03 CONVERT GAUSSIAN OUTPUT TO MOLFILE

The next step is to produce a hybrid molfile v2000/3000 file. The program babble we built in the first step can convert the gaussian output to a molfile that has the v2000 structure but uses the v3000 bond types. In particular the molfile_2_params_polymer.py script needs to have aromatic bonds types (type 4) for aromatic rings instead of alternating single and double bonds (Kekule structure). 

The command bellow will convert the gaussian output to molfile format for the example...


```
> openbabel/install/bin/babel -i g03 stage_03_gaussian_output/ornithine.log -o mol stage_04_molfile/ornithine.mol
```

```
<<<<< ornithine.mol start >>>>>
 OpenBabel01101117083D

 31 30  0  0  0  0  0  0  0  0999 V2000
    0.4866    2.2296   -0.2354 C   0  0  0  0  0
    1.2081    1.9315   -1.1549 O   0  0  0  0  0
    0.4369    3.6268    0.3352 C   0  0  0  0  0
   -0.3187    4.1933   -0.1996 H   0  0  0  0  0
    0.1852    3.6362    1.3888 H   0  0  0  0  0
    1.3920    4.1075    0.1799 H   0  0  0  0  0
   -2.6920   -1.1572   -0.6650 N   0  0  0  0  0
   -4.0550   -1.5847   -0.3922 C   0  0  0  0  0
   -2.3500   -1.2554   -1.5934 H   0  0  0  0  0
   -4.1238   -1.9440    0.6231 H   0  0  0  0  0
   -4.7601   -0.7727   -0.5220 H   0  0  0  0  0
   -4.3066   -2.3874   -1.0709 H   0  0  0  0  0
   -0.3415    1.3293    0.3431 N   0  0  0  0  0
   -0.6023    0.0286   -0.2270 C   0  0  0  0  0
   -2.0235   -0.3660    0.1877 C   0  0  0  0  0
   -2.4547   -0.0133    1.2526 O   0  0  0  0  0
    0.3558   -1.0672    0.2829 C   0  0  0  0  0
    1.8184   -0.8036   -0.0848 C   0  0  0  0  0
    4.1559   -1.6146    0.0004 N   0  0  0  0  0
    2.7190   -1.9434    0.3678 C   0  0  0  0  0
   -0.9580    1.6065    1.0765 H   0  0  0  0  0
   -0.5197    0.1005   -1.3044 H   0  0  0  0  0
    0.2486   -1.1369    1.3608 H   0  0  0  0  0
    0.0400   -2.0197   -0.1352 H   0  0  0  0  0
    1.9103   -0.6656   -1.1566 H   0  0  0  0  0
    4.2542   -1.4797   -0.9970 H   0  0  0  0  0
    2.1487    0.1188    0.3755 H   0  0  0  0  0
    4.7999   -2.3436    0.2758 H   0  0  0  0  0
    2.4931   -2.8811   -0.1191 H   0  0  0  0  0
    2.7087   -2.0883    1.4385 H   0  0  0  0  0
    4.4568   -0.7571    0.4437 H   0  0  0  0  0
  2  1  2  0  0  0
  1  3  1  0  0  0
  1 13  1  0  0  0
  4  3  1  0  0  0
  6  3  1  0  0  0
  3  5  1  0  0  0
  9  7  1  0  0  0
  7  8  1  0  0  0
  7 15  1  0  0  0
 12  8  1  0  0  0
 11  8  1  0  0  0
  8 10  1  0  0  0
 14 13  1  0  0  0
 13 21  1  0  0  0
 22 14  1  0  0  0
 14 15  1  0  0  0
 14 17  1  0  0  0
 15 16  2  0  0  0
 24 17  1  0  0  0
 18 17  1  0  0  0
 17 23  1  0  0  0
 25 18  1  0  0  0
 18 20  1  0  0  0
 18 27  1  0  0  0
 26 19  1  0  0  0
 19 28  1  0  0  0
 19 20  1  0  0  0
 19 31  1  0  0  0
 29 20  1  0  0  0
 20 30  1  0  0  0
M  END
<<<<< ornithine.mol end >>>>>
```

-----------------------------------
   STEP 04 MODIFYING THE MOLFILES
-----------------------------------
The molfile2params_polymer.py script requires some additional data to be added to the end of the molfile. This data is specified at the end of the file after the bond information. It is a list of variable names and then a list of values. The variable are described bellow. 

ROOT: Single numerical value. Atom number (according to the order in the molfile). Where the atom tree is rooted for this residue type. Should be the nitrogen of the central residue.
POLY_N_BB, POLY_CA_BB, POLY_C_BB, POLY_CO_BB: Atom number (according to the order in the molfile). The backbone nitrogen, alpha-carbon, carbonyl-carbon and carbonyl-oxygen. These get special rosetta atom types and so are listed here special.
POLY_IGNORE: List of atom numbers (according to the order in the molfile). These are the atoms for the capping groups with the exception of the upper and lower connect atoms. They will not be listed in the atoms and bonds in the params file but are used in determining the atom types.
POLY_UPPER, POLY_LOWER: Atom number (according to the order in the molfile). These are the atoms in the capping groups that connect to the residue. They will not be listed in the atoms and bonds in the params file but are used in determining the atom types and they are listed in the internal coordinate section.
POLY_CHG: Single numerical value. Overall charge on the residue. 
POLY_PROPERTIES: List of alpha-numerical values. These get used by Rosetta at various places in the program. You can say something like 

```
if ( pose.residue(10).type().is_protein() ) { // do something }.
END: The end of the file.
```


Note: There are 2 spaces between the "M" and the variable name.
Note: If you have to make multiple residue types keeping the atoms in the same order makes assigning all these numbers easier. 

Additional info for peptide ornithine...

```
M  ROOT 13
M  POLY_N_BB 13
M  POLY_CA_BB 14
M  POLY_C_BB 15
M  POLY_O_BB 16
M  POLY_IGNORE 2 3 4 5 6 8 9 10 11 12
M  POLY_UPPER 7
M  POLY_LOWER 1
M  POLY_CHG 1
M  POLY_PROPERTIES PROTEIN POLAR CHARGED
M  END
```
----------
###   STEP 05 RUNNING MOLFILE2PARAMS.PY

Finally, the molfile2params_polymer.py script will convert the modified molfile to a params file. The commands bellow work for the example files and produce a ResidueType parameter file for ornithine with the three letter code C40. 

```
$python scripts/molfile_to_params_polymer.py --clobber --polymer --no-pdb --name C40 -k ornithine.kin stage_05_modified_molfile/ornithine.mol
```


There may need additional tweaking that needs to happen to the params files to make them work correctly. Compare these to the ones in the database for reference. In particular it is important to double check that the chi angles and the atoms that comprise each of them are specified correctly. Additionally it is important to check that the dependancies in the internal coordinates are correctly setup with respect the the atoms specified for the chi angles. For example, in ornithine the chi1 is defined as

`CHI 1  N    CA   CB   CG` 

and the internal coordinates for it and it hydrogens are

```
ICOOR_INTERNAL    CG   -62.608549   66.937720    1.530976   CB    CA    N  
ICOOR_INTERNAL   1HG   121.953676   70.110328    1.084548   CG    CB    CD 
ICOOR_INTERNAL   2HG   116.869954   70.361318    1.082495   CG    CB   1HG 
```

Since the gamma hydrogens are defined relative to the CG they will move when the chi1 is rotated.

Pro Tip: The molfile2params.py script can produce a kinemage file. This is handy as it lets you check how rosetta will build the atom tree, and the atom type assignments, and other parameters. You can open kinemage files using the KiNG program from the Richardson lab at Duke (http://kinemage.biochem.duke.edu/software/king.php).

#### Use the new params file:
The generated params file can be used in protocols using the `-extra_res_fa` or `-extra_res_cen` options.

-

## Part II - CALCULATING EXPLICIT UNFOLDED STATE ENERGIES FOR NONCANONICAL AMINO ACIDS

KEYWORDS: NONCANONICALS ANALYSIS
Calculating the explicit unfolded state energies is the third of three steps toward being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles and HowToMakeRotamerLibraries. 

The explicit unfolded state energies of an amino acid represent the energy of an amino acid in the unfolded state of a protein and is used to replace the reference energies in Rosetta. The UnfoldedStateEnergyCalculator uses a fragment based method to calculate the average unfolded state energies for each ResidueType. The protocols works on a large set of protein structures that are split in to randomly generated fragments. The central residue of each fragment is mutated to the residue of interest. The fragment is repacked. The unweighted energy for each energy method in the scoring function is recorded for the central residue. After the energies for all fragment central residues are collected, a boltzmann-weighted-average average energy is calculated for each term. 

Calculation explicit unfolded state energies for a NCAA requires three steps:
 - Obtaining a set of input pdbs
 - Running the UnfoldedStateEnergyCalculator protocol on the set of pdbs
 - Modifying the unfolded state energies file in the database

-

Some preparation first:

```
$> ln -s  HowToMakeExpliciteUnfoldedStateEnergies/inputs .
$> ln -s  HowToMakeExpliciteUnfoldedStateEnergies/scripts .
```

## How To Make Explicite Unfolded State Energies


###   STEP 01 OBTAINING A SET OF INPUT PDBS


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

`ls inputs | grep '\.pdb' | while read line; do echo "inputs/$line" >> inputs/list; done` 


**Citation: G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling server. Bioinformatics, 19:1589-1591, 2003. 


###   STEP 02 RUNNING THE UNFOLDEDSTATEENERGYCALCULATOR PROTOCOL


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
$> $ROSETTA3/UnfoldedStateEnergyCalculator.linuxclangrelease -ignore_unrecognized_res -extrachi_cutoff 0 -l inputs/list -residue_name C40 -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_c40
```

NOTE: The extension on your executable my be different.

The run will take between 30-60 seconds per pdb file.

The log file contains lots of useful information. It contains the unweighted energies for each of the energy methods for each of the individual fragments. At the end it will print the average unweighted energies for each ResidueType as well as the Boltzmann weighted average unweighted energies. Boltzmann weighted average unweighted energies are used because some backbones just can't tolerate a mutation to a particular ResidueType and there are extremely high repulsive energies for some fragments that skew the average value. Using the Boltzmann weighting removes the higher energy outliers in a more elegant fashion than a hard energy cutoff.

--------------------------------------------
###   STEP 03 MODIFY THE UNFOLDED STATE ENERGIES FILE

Once the UnfoldedStateEnergyCalculator has finished running the Boltzmann weighted average unweighted energies need to be added to the database. The line you want is the "BOLZMANN UNFOLDED ENERGIES". These are the Boltzmann weighted average unfolded energies for each energy method. The file you need to modify is unfolded_state_residue_energies_mm_std.

Using the ornithine line as an example, the line form the log file is... 

BOLZMANN UNFOLDED ENERGIES:  fa_atr:    -2.462 fa_rep:     1.545 fa_sol:     1.166 mm_lj_intra_rep:     1.933 mm_lj_intra_atr:    -1.997 mm_twist:     2.733 pro_close:     0.009 hbond_sr_bb:    -0.006 hbond_lr_bb:     0.000 hbond_bb_sc:    -0.001 hbond_sc:     0.000 dslf_ss_dst:     0.000 dslf_cs_ang:     0.000 dslf_ss_dih:     0.000 dslf_ca_dih:     0.000

We could add the following to the unfolded\_state\_residue\_energies\_mm_std file in the database using the command bellow.

```
$ echo "C40 -2.462 1.545 1.166 1.933 -1.997 2.733 0.009 -0.006 0.000 -0.001  0.000" >> ROSETTA/main/database/scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std 
```

The ResidueType can now be used in almost any Rosetta protocol that is compatible with the MM_STD scoring function.



## Part IV    INSTRUCTIONS FOR RUNNING THE DOUGSDOCKDESIGNMINIMIZE PROTOCOL


KEYWORDS: DOCKING INTERFACES STRUCTURE_PREDICTION

The DougsDockDesignMinimize (DDDM) protocol was used in the accompanying manuscript to redesign the protein/peptide interface of Calpain and a fragment of its inhibitory peptide calpastatin. The protocol was written for this specific protein/peptide interaction and modifications to the code will be necessary to run the protocol on a different system. A modified form was used as the example protocol in the advanced section of the Rosetta 3.0 release manual.

Minor modifications to the protocol have been made from the version used to produce the designs in the accompanying publication. The interfaces of the Rosetta libraries have changed since the initial implementation of the protocol and these modifications were necessary to allow the protocol to work with the current release. 

The protocol has two loops referred to as the inner and outer loops. The outer loop controls the number of structures generated, outputting structures only if certain filters are met. The inner loop iterates between two phases: a perturbation phase and a design phase. During the perturbation phase three types of perturbations are used. Two types of perturbations are applied to the peptide: rotational and translational rigid body permutation of the peptide in the binding pocket of the protein, and small and shear perturbations of the backbone dihedral angles in all residues of the peptide. The third type of perturbation is preformed on residues 1-45 of the N-terminus which comprise the peptide binding site and surrounding residues of the protein are perturbed using small and shear perturbations of the backbone dihedral angles. One of the three perturbations is randomly chosen to be used each time the function in called and is followed by a coarse side chain optimization (RotamerTrials). During the design phase of the protocol the side chains positions are optimized using a more intensive side chain optimization routine (PackRotamers) followed by minimization of the backbone and side chain dihedrals angles as well as the jump between the peptide and protein. The number of perturbations per design as well as the magnitude of perturbations are controlled using command line options described bellow.

There are three basic steps to running the protocol as it was used in the accompanying manuscript.
 - generating the input resfiles and folders
 - running the DougsDocDesignMinimize protocol for each mutation at each position 
 - running the analysis scripts

Example output is provided for a small version of the full run.

------------------------------------------
   GENERATING INPUT RESFILES AND FOLDERS
------------------------------------------
The Protein Databank code for the Calpain/Calpastatin structure used in the designs is 1NX1. The files contain two copies of calpain (chains A and B) and two copies of the calpastatin inhibitory peptide (chains C and D). To reduce the size of the system, designs were done on chain A and chain C only because chains A and C have lower B-factors than chains B and D. The calpain/calpastatin interface is distal to the calpain/calpain interface and the resides that make up the calpain/calpain interface were held fix during the protocol. Additionally the calcium atoms and water molecules were removed. Rosetta doesn't work with with water molecules and the calcium ions are not near the protein peptide interface. This input pdb was repacked and both the side chain and backbone dihedrals minimized using a modified version of the fixbb app. The input PDB file is called 1NX1_clean_repack_min_all.pdb

Each of the NCAAs added in the accompanying publication was tried at each position in the peptide. To do this and to keep all of the output pdb organized a script is provided that creates folders and generates resfiles based on templates for each sequence. To generate the resfiles and folders run the following commands...

```
$> ln -s HowToRunDougsDockDesignMinimizeProtocol/scripts .
$> ln -s HowToRunDougsDockDesignMinimizeProtocol/run_dir/* .
$> scripts/make_folders_resfiles.tcsh
```

NOTE: To save space, the script has been modified to only produce resfiles and folders for position 610 in the peptide and residue type MPA (4-methyl-phenylalanine). To modify the script to produce folders and resfiles for each position simply uncomment the lines that read "#foreach i ( 601 602 603 604 605 606 607 608 609 610 611 )" and "#foreach j ( ABA APA HLU ..... C92 C93 C94 )" comment out the line that reads "foreach i ( 610 )" and "foreach j ( MPA )".

Using the templates in the run_dir the make_folders_resfiles.tcsh script makes folders and resfiles to preform all of the DDDMI runs.



###   RUNNING DOUGSDOCKDESIGNMINIMIZEINTERFACE PROTOCOL


To run the protocol modifications need to made to files in the database. In the file rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt all of the paths to the residue type parameter files under the L-NCAA heading need to be uncommented by removing the "#" from the front of the line. Additionally the rotamer libraries for the NCAA are not provided in the default Rosetta database because they are more than 400MB. The rotamer libraries for the NCAAs added in the accompanying publication are provided as supplemental information. 

NOTE: Turning on all of the additional residue types dramatically increases the number of residue types and the memory footprint of Rosetta. The memory foot print can be reduced by commenting out unnecessary patches in the rosetta_database/chemical/residue_type_sets/fa_standard/patches.txt file. For the DougsDockDesignMinimizeProtocol all but the NtermProteinFull.txt and CtermProteinFull.txt can be safely commented out by placing "#" symbols at the beginning of each line of the patches.txt except for the lines that say "NtermProteinFull.txt" and "CtermProteinFull.txt".

To run the protocol as in the accompanying publication preform the following commands starting at the HowToRunDougsDockDesignMinimizeProtocol directory.

```
$> ln -s pos_610_MPA/* .
$> $ROSETTA/bin/doug_dock_design_min_mod2_cal_cal.linuxgccrelease -database /PATH/TO/rosetta_database -s inputs/1NX1_clean_repack_min_all.pdb -resfile resfile_pos_603_MPA -nstruct 1 -inner_num 45 -pert_num 25 -ia_ener 100 -use_input_sc -pdb_gz
```

>**change nstruct to 255!**

The above command generates 255 structures and will take approximately 5 minutes per structure depending on your hardware.

NOTE: The extension of your executable maybe different than the above. Also in the publication the nstruct command line option was 255. However to save space for the protocol capture an nstruct of 10 was used.

A script is provided that will preform the above command for each folder created by the make_folders_resfiles.tcsh script.

```
$> for i in pos_*; do cd $i $ROSETTA3/doug_dock_design_min_mod2_cal_cal.linuxgccrelease ../inputs/ \  
1NX1_clean_repack_min_all.pdb -resfile ../resfile_$i \  
-nstruct 1 -inner_num 45 -pert_num 25 -ia_ener 100 \  
-use_input_sc -pdb_gz >& $i.log; \ 
 echo "Finished $i"; cd ../; done
          
```

NOTE: You will need to set the path to your database and executable in the run_script.bash.

There are a number command line options to control the running of the application.

pert_mc_temp: The MC temperature to use for the perturbation phase of the DDDM protocol, defaults to 0.8 kT.
pert_dock_rot_mag: The rotation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.5. 
pert_dock_trans_mag: The translation magnitude for the ridged body perturbation in the perturbation phase of the DDDM protocol, defaults to 0.25. 
pert_pep_small_temp: The temperature for the internal MC object in the small mover of the peptide perturbations, defaults to 0.8 kT. 
pert_pep_shear_temp: The temperature for the internal MC object in the shear mover of the peptide perturbations, defaults to 0.8 kT. 
pert_ter_small_temp: The temperature for the internal MC object in the small mover of the termini perturbations, defaults to 0.8 kT. 
pert_ter_shear_temp: The temperature for the internal MC object in the shear mover of the termini perturbations, defaults to 0.8 kT. 
pert_pep_small_H: The maximum angle of perturbation for helical secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_small_L: The maximum angle of perturbation for loop secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_small_E: The maximum angle of perturbation for strand secondary structure for the peptide small mover, defaults to 1 degree.
pert_pep_shear_H: The maximum angle of perturbation for helical secondary structure for the peptide shear mover, defaults to 1 degree
pert_pep_shear_L: The maximum angle of perturbation for loop secondary structure for the peptide shear mover, defaults to 1 degree.  
pert_pep_shear_E: The maximum angle of perturbation for strand secondary structure for the peptide shear mover, defaults to 1 degree.
pert_ter_small_H: The maximum angle of perturbation for helical secondary structure for the termini small mover, defaults to 1 degree
pert_ter_small_L: The maximum angle of perturbation for loop secondary structure for the termini small mover, defaults to 1 degree.  
pert_ter_small_E: The maximum angle of perturbation for strand secondary structure for the termini small mover, defaults to 1 degree.
pert_ter_shear_H: The maximum angle of perturbation for helical secondary structure for the termini shear mover, defaults to 1 degree
pert_ter_shear_L: The maximum angle of perturbation for loop secondary structure for the termini shear mover, defaults to 1 degree.  
pert_ter_shear_E: The maximum angle of perturbation for strand secondary structure for the termini shear mover, defaults to 1 degree.
pert_pep_num_rep: Number of small and shear iterations for the peptide, defaults to 100. 
pert_ter_num_rep: Number of small and shear iterations for the terminus, defaults to 100. 
pert_num: Number of iterations of perturbation loop per design, defaults to 100.
inner_num: Number of iterations of the inner loop, defaults to 100. 
ia_ener: Upper energy limit for final design/interface analysis checkpoint, defaults to 0.0. 
desn_mc_temp: The temperature to use for the design/minimization phase of the DDDM protocol, defaults to 0.8 kT.

For the most part the defaults should suffice and are what was used in the paper. The "ia_ener" command line option is dependent on the complex that the protocol is run on. Setting it unreasonably low will cause the protocol to run forever and never output any structures. Setting it to 100 above allows all but the most aggressions structures to be output. The nstruct command line option controls the outer loop. 

-----------------------------
   ANALYZING THE RESULTS
-----------------------------

Each of the output pdb files contains information about the protein peptide complex that can used to evaluate the designs.  For example at the end of the 1NX1_clean_repack_min_all_0004.pdb is shown bellow. For each filter the value is calculated for the protein and peptide together (COMPLEX) and separated by 1000 angstroms (SEPARATE) and the difference between the two (DIFF). ENERGY is the Rosetta energy. The ENERGY_COMPLEX is the primary determinant to how good a design is and the ENERGY_DIFF can give an estimate for the binding energy. SASA is the solvent accessible surface area. SASA_DIFF is indicative of sequences that make a more protein-peptide contacts and can be used for screening designs for example placing a very large side chain at a constrained interface position can cause the peptide to be pushed out of the binding pocket which would be reflected in a smaller magnitude SASA_DIFF. HB_ENER is the hydrogen bonding component of the Rosetta energy. Larger HB_ENER_DIFF values indicate that the design is making more or better hydrogen bonds across the protein peptide interface. PACK is the RosettaHoles score and is a measurement of how well the protein is packed. A PACK_COMPLEX that is larger than the PACK_SEPARATE is favorable and suggests that the complex is better packed than the protein alone. Additionally the RosettaHoles score penalizes holes that cannot be occupied by solvent so larger PACK_DIFF score indicate that the designed peptide is capable of filling cavities in the protein that are inaccessible to solvent.

```
ENERGY_COMPLEX:	   -63.3501
ENERGY_SEPERATE:   -48.599
ENERGY_DIFF:	   -14.7512
SASA_COMPLEX:	   11774.1
SASA_SEPERATE:	   13092.3
SASA_DIFF:	   -1318.15
HB_ENER_COMPLEX:   -146.728
HB_ENER_SEPERATE:  -145.037
HB_ENER_DIFF:	   -1.69085
PACK_COMPLEX:	   0.515277
PACK_SEPERATE:	   0.466995
PACK_DIFF:	   0.0482824
```

A script is provided that pulls the information out of a set of pdb files and sorts it based on the ENERGY_COMPLEX metric. 

```
$ cd pos_610_MPA
$ ../scripts/get_interface_data.tcsh
```

