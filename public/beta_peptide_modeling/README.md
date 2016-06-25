Beta-3-peptide design
=====================

KEYWORDS: PEPTIDES DESIGN

This demo illustrates a protocol to redesign and model beta-3-peptides in Rosetta.
It was written in July 2012 by Fang-Chieh Chou (fcchou at stanford dot edu).
Detailed application of the method is described in the following paper:

* Molski, M.A., Goodman, J.L., Chou, F.C., Baker, D., Das, R., and Schepartz, A. (2012). Remodeling a beta-peptide bundle. Chemical Science.

Running the Demo
----------------

The example input files are in rosetta_input.

The example command lines are given in rosetta_input/cmdline.
Cd into rosetta_input/ and type `sh cmdline` to run the commands.
After the job finished successfully, you should see the output pdb files in the corresponding folders and log files in the current folder.
Sample outputs are given in the example_output folder for comparison.

Command-lines
-------------

1. Redesigning

```bash
$> $ROSETTA3/bin/beta_peptide_modeling.default.linuxgccrelease -force_field beta_peptide -native rosetta_inputs/redesign/acdy_LLLL_LLLL.pdb -algorithm redesign -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 -repack_res 2 5 8 11 14 17 20 23 -n_repeat 4 -repeat_size 24
```

   `-native` gives the starting model.

   `-repack_res` specify the residue being redesigned/repacked.

   `-n_repeat` and `-repeat_size` set up the symmetry based repacking. In this example, each residue n listed in `repack_res` and n+24, n+48, n+72 (4 copies total) are repacked/redesigned. Also each set of 4 residues are enforced to be the same residue type and have the same rotameric conformation.

2. Repacking

```bash
$> $ROSETTA3/bin/beta_peptide_modeling.default.linuxgccrelease -database ../../../../rosetta_database -force_field beta_peptide_soft_rep_design -native rosetta_inputs/repack_and_minimize/acdy_LFFL_LFFL.pdb -algorithm repack -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 -repack_res 2 5 8 11 14 17 20 23 -n_repeat 4 -repeat_size 24
```
   Similar to above, but only repack the original residues.

3. Minimization

```bash
$> $ROSETTA3/bin/beta_peptide_modeling.default.linuxgccrelease -database  ../../../../rosetta_database -force_field beta_peptide -native rosetta_inputs/repack_and_minimize/acdy_LFFL_LFFL_repack.pdb -algorithm minimize -ex1 -ex2 -packing::pack_missing_sidechains false -packing::extrachi_cutoff 0 
```

   Minimize a given pdb, return the final model and Rosetta energies.

Extra flags and score files
---------------------------

`-score::no_smooth_etables` true Use the old Rosetta attractive/repulsive score without smoothing.

`-no_symmetry` true Disable the symmetry-based redesigning.

`-force_field beta_peptide_mm` / `-force_field beta_peptide_soft_rep_mm` Use the molecular mechanics-based Rosetta scoring files. Replace the corresponding `-force_field` flags with the MM ones in the above command lines. Also, extra flags `-apply_dihedral_cst false` should be used to turn off dihedral constraint when MM potential is used.

