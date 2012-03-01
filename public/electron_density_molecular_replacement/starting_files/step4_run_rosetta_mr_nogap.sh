#!/bin/sh

mkdir relaxes
cd relaxes

$ROSETTA3_SRC/bin/mr_protocols.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-MR::mode cm \
	-in::file::extended_pose 1 \
	-in::file::fasta ../inputs/1crb.fasta \
	-in::file::alignment ../templates/2qo4.ali \
	-in::file::template_pdb ../phaser/2qo4_mr.PHASER.1.pdb \
	-edensity:mapreso 3.0 \
	-edensity:grid_spacing 1.5 \
	-edensity:mapfile ../phaser/2qo4_mr.PHASER.1_2mFo-DFc.ccp4 \
	-edensity::sliding_window_wt 1.0 \
	-edensity::sliding_window 5 \
	-relax::default_repeats 2 \
	-relax::jump_move true \
	-cm::aln_format grishin \
	-MR::max_gaplength_to_model 0 \
	-nstruct 1 \
	-ignore_unrecognized_res -overwrite -out:prefix nogaps_
