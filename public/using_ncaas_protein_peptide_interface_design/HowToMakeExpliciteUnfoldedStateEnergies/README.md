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

Since the UnfoldedStateEnergyCalculator protocol uses fragments from protein structures, we need a set of high quality structures to work with. Through their PISCES server, the Dunbrack laboratory maintains lists of structures in the Protein Data Bank organized based on xray resolution, precent sequence similarity, and r-factors**. These lists are a convenient way to get a set of high quality structures. In this example we will use a list culled on May 20, 2011. It contains 1801 pdb files that have an an xray resolution of at least 1.6 angstroms, less than 20% sequence identity, and r-factors of less than 0.25. To get the pdbs simply use a supplied script to download the pdbs from the Protein Data Bank ftp servers. 

$> cd inputs
$> ../scripts/get_pdbs.bash cullpdb_pc20_res1.6_R0.25_d110520_chains1859

There should be 1801 gzipped pdb files and a text file containing a list of them called cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned in the inputs directory. Rosetta will sometimes fail to correctly read in particular pdbs files. The cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned file is a list of the pdbs which have been screened to be read successfully by Rosetta. 

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

$> cd outputs
$> $ROSETTA3/UnfoldedStateEnergyCalculator.linuxclangrelease -ignore_unrecognized_res -extrachi_cutoff 0 -l ../inputs/cullpdb_pc20_res1.6_R0.25_d110520_chains1859_list_pruned -residue_name C40 -mute all -unmute devel.UnfoldedStateEnergyCalculator -unmute protocols.jd2.PDBJobInputer -no_optH true -detect_disulf false >& ufsec_log_c40 &

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
