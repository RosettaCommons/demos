# Instructions for Calculating Explicit Unfolded State Energies for Noncanonical Amino Acids

KEYWORDS: ANALYSIS GENERAL

This README was written by P. Douglas Renfrew (renfrew@nyu.edu). The file was updated by Parisa Hosseinzadeh (parisah@uw.edu).

This demo illustrates the UnfoldedStateEnergyCalculator protocol. It was originally written as the 3rd (of 4) parts of a "protocol capture" accompanying the paper "Using Noncanonical Amino Acids in Computational Protein-Peptide Interface Design" by P. Douglas Renfrew, Eun Jung Choi, Brian Kuhlman (2011) PLoS One, from the RosettaCon2010 special PLoS ONE special issue. The four separate protocol captures that describe the four different ways Rosetta was used in the publication: creating Rosetta ResidueType parameter files for NCAAs, creating backbone dependent rotamer libraries for NCAAs, calculating explicit unfolded state reference energies for NCAAs, and the running of the DougsDockDesignMinimizeInterface protocol using NCAAs. The protocol captures for creating ResidueType parameter files, rotamer libraries and explicit unfolded state energies describe the process used in the publication but are written from the standpoint of a researcher looking to add an additional NCAA to Rosetta. 

Calculating the explicit unfolded state energies is the third of three steps toward being able to use a noncanonical amino acid (NCAA) in Rosetta. To add a new NCAA or to better understand how the NCAAs in the related publication were added one should have already completed or understand the steps in HowToMakeResidueTypeParamFiles and HowToMakeRotamerLibraries. 

The explicit unfolded state energies of an amino acid represent the energy of an amino acid in the unfolded state of a protein and is used to replace the reference energies in Rosetta. The UnfoldedStateEnergyCalculator uses a fragment based method to calculate the average unfolded state energies for each ResidueType. The protocols works on a large set of protein structures that are split in to randomly generated fragments. The central residue of each fragment is mutated to the residue of interest. The fragment is repacked. The unweighted energy for each energy method in the scoring function is recorded for the central residue. After the energies for all fragment central residues are collected, a boltzmann-weighted-average average energy is calculated for each term. 

Calculation explicit unfolded state energies for a NCAA requires three steps:
 - Obtaining a set of input pdbs
 - Running the UnfoldedStateEnergyCalculator protocol on the set of pdbs
 - Modifying the unfolded state energies file in the database

# Step 1: Obtaining a Set of Input PDBs

Since the UnfoldedStateEnergyCalculator protocol uses fragments from protein structures, we need a set of high quality structures to work with. Through their PISCES server, the Dunbrack laboratory maintains lists of structures in the Protein Data Bank organized based on xray resolution, precent sequence similarity, and r-factors\*\*. These lists are a convenient way to get a set of high quality structures. In this example we will use a list culled on May 20, 2011. It contains 1801 pdb files that have an xray resolution of at least 1.6 angstroms, less than 20% sequence identity, and r-factors of less than 0.25. To get the pdbs simply use a supplied script to download the pdbs from the Protein Data Bank ftp servers. 

```
$ cd inputs
$ ../scripts/get_pdbs.bash cullpdb_pc20_res1.6_R0.25_d110520_chains1859
```
There should be 1801 gzipped pdb files and a text file containing a list of them called cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned in the inputs directory. Rosetta will sometimes fail to correctly read in particular pdbs files. The cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned file is a list of the pdbs which have been screened to be read successfully by Rosetta. 

\*\*Citation: G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling server. Bioinformatics, 19:1589-1591, 2003. 

# Step 2: Running the UnfoldedStateEnergyCaclulator Protocol

The UnfoldedStateEnergyCalculator is relatively easy to run. The command line options are described below:

- frag_size: single integer value, sets the number of residues in each fragment, should be an odd number and has a default of 5 which is what was used in the accompanying publication
- residue_name: string value, sets the three letter code of the residue type which the central residue will be mutated to
- repack_fragments: boolean value, controls if the fragments will be repacked before scoring and defaults to true
- native_sequence: boolean value, controls if the central residue will be mutated before scoring and defaults to false

Additionally it is strongly recommended to add the following flags as they will make Rosetta handle more pdb files and improves runtime by disabling default features that will be negated by the fragmenting and prepacking:

- ignore_unrecognized_res: causes Rosetta to ignore unrecognized residue types and 
- ex1 and ex2 and extrachi_cutoff 0: force rosetta to use additional rotamer during the fragment repacking
- mute all and unmute devel.UnfoldedStateEnergyCalculator and unmute protocols.jd2.PDBJobInputer: reduces the size of the log file significantly by turning off unnecessary output
- no_optH true: turns off the hydrogen optimization done when the protein is first read in 
- detect_disulf false: turns off disulfide detection

Continuing the ornithine example we have used in the two previous protocol captures, to calculate the unfolded state energies one would run the following command.

```
$ cd outputs
$ PATH/TO/bin/UnfoldedStateEnergyCalculator.linuxgccrelease -ignore_unrecognized_res -ex1 -ex2 -extrachi_cutoff 0 -l ../inputs/cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned -residue_name C40 -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_c40 &
```

**Note:** The extension on your executable my be different.
**Note:** if you don't have the C40 .params and .rotlib from previous steps, this code will fail. You can either test it with your ncAA files or use NVL instead of C40 in the code.

The run will take between 30-60 seconds per pdb file.

The log file contains lots of useful information. It contains the unweighted energies for each of the energy methods for each of the individual fragments. At the end it will print the average unweighted energies for each ResidueType as well as the Boltzmann weighted average unweighted energies. Boltzmann weighted average unweighted energies are used because some backbones just can't tolerate a mutation to a particular ResidueType and there are extremely high repulsive energies for some fragments that skew the average value. Using the Boltzmann weighting removes the higher energy outliers in a more elegant fashion than a hard energy cutoff.

# Step 3: Modify the Unfolded State Energies File

Once the UnfoldedStateEnergyCalculator has finished running the Boltzmann weighted average unweighted energies need to be added to the database. The line you want is the "BOLZMANN UNFOLDED ENERGIES". These are the Boltzmann weighted average unfolded energies for each energy method. The file you need to modify is unfolded_state_residue_energies_mm_std.

Using the ornithine line as an example, the line form the log file is... 

```
BOLZMANN UNFOLDED ENERGIES:  fa_atr:    -2.462 fa_rep:     1.545 fa_sol:     1.166 mm_lj_intra_rep:     1.933 mm_lj_intra_atr:    -1.997 mm_twist:     2.733 pro_close:     0.009 hbond_sr_bb:    -0.006 hbond_lr_bb:     0.000 hbond_bb_sc:    -0.001 hbond_sc:     0.000 dslf_ss_dst:     0.000 dslf_cs_ang:     0.000 dslf_ss_dih:     0.000 dslf_ca_dih:     0.000
```

We could add the following to the unfolded_state_residue_energies_mm_std file in the database using the command bellow.

```
$ echo "C40 -2.462 1.545 1.166 1.933 -1.997 2.733 0.009 -0.006 0.000 -0.001  0.000" >> minirosetta_database/scoring/score_functions/unfolded/unfolded_state_residue_energies_mm_std 
```

The ResidueType can now be used in almost any Rosetta protocol that is compatible with the MM_STD scoring function.

# Shorter runs

If you want to test how things work, you can test a shorter version of the codes mentioned.

go to the short run files:
```
$ cd short_run_files
```
Then obtain the PDBs running:
```
$ ../scripts/get_pdbs.bash small_list
```
Alternatively, you can copy this file to your directory:
```
$> cp short_run_files/small_list_pruned .
```

Next step is to run the calculator using this command:(`$ROSETTA3`=path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/UnfoldedStateEnergyCalculator.default.linuxgccrelease -ignore_unrecognized_res -ex1 -ex2 -extrachi_cutoff 0 -l small_list_pruned -residue_name NVL -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_NVL_short
```
The expected outputs are stored in the output_short directory inside the short_run_files.
