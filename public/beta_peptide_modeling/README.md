####################
 Beta-3-peptide designing 
####################

This README was written in Jul. 2012, by Fang-Chieh Chou (fcchou@stanford.edu).

This demo illustrates a protocol to redesign and model beta-3-peptides in Rosetta.

Detailed application of the method is described in the following paper:

Molski, M.A., Goodman, J.L., Chou, F.C., Baker, D., Das, R., and Schepartz, A. (2012). Remodeling a beta-peptide bundle. Chemical Science.

The example input files are in rosetta_input.

The beta-3-peptide residues is less commonly used so is not turned on by default. Therefore the following setup steps are required prior to running this demo.

Before running the command lines, edit the file rosetta_database/chemical/residue_type_sets/fa_standard/residue_types.txt, find the section "Beta-peptide Types" and uncommented all the items listed below (residue_types/beta-peptide/B3A.params etc.).

The example command lines are given in rosetta_input/cmdline. Cd into rosetta_input/ and type "sh cmdline" to run the commands.

After the job finished sucessfully, you should see the output pdb files in the corresponding folders and log files in the current folder. Sample outputs are given in the example_output folder for comparsion.

Detailed explanation of the command lines are given below.

1. Redesigning
beta_peptide_modeling.linuxgccrelease -database  ../../../../rosetta_database -force_field beta_peptide -native redesign/acdy_LLLL_LLLL.pdb -algorithm redesign -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 -repack_res 2 5 8 11 14 17 20 23 -n_repeat 4 -repeat_size 24

-native gives the starting model.
-repack_res specify the residue being redesigned/repacked.
-n_repeat and -repeat_size set up the symmetry based repacking. In this example, each residue n listed in "repack_res" and n+24, n+48, n+72 (4 copies total) are repacked/redesigned. Also each set of 4 residues are enforced to be the same residue type and have the same rotameric conformation.

2. Repacking
beta_peptide_modeling.linuxgccrelease -database ../../../../rosetta_database -force_field beta_peptide_soft_rep_design -native repack_and_minimize/acdy_LFFL_LFFL.pdb -algorithm repack -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 -repack_res 2 5 8 11 14 17 20 23 -n_repeat 4 -repeat_size 24
Similar to above, but only repack the original residues.

3. Minimization
beta_peptide_modeling.linuxgccrelease -database  ../../../../rosetta_database -force_field beta_peptide -native repack_and_minimize/acdy_LFFL_LFFL_repack.pdb -algorithm minimize -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 
Minimize a given pdb, return the final model and Rosetta energies.

Extra flags / scoring files

-score::no_smooth_etables true Use the old Rosetta attractive/repulsive score without smoothing.

-no_symmetry true Disable the symmetry-based redesigning.

-force_field beta_peptide_mm / -force_field beta_peptide_soft_rep_mm Use the molecular mechanics-based Rosetta scoring files. Replace the corresponding “-force_field” flags with the MM ones in the above command lines. Also, extra flags “-apply_dihedral_cst false” should be used to turn off dihedral constraint when MM potential is used.

