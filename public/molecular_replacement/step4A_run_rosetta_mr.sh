#!/bin/sh

mkdir -p outputs
cd outputs

../../../../main/source/bin/mr_protocols.default.linuxclangrelease\
	-in::file::fasta ../inputs/1crb.fasta \
	-in::file::alignment ../templates/2qo4.ali \
	-in::file::template_pdb ../phaser/2qo4_mr.PHASER.1.pdb \
 	-loops::frag_files ../inputs/xx1crb_09_05.200_v1_3.gz ../inputs/xx1crb_03_05.200_v1_3.gz none \
	-edensity:mapreso 3.0 \
	-edensity:grid_spacing 1.5 \
	-edensity:mapfile ../phaser/2qo4_mr.PHASER.1_2mFo-DFc.ccp4 \
	-MR::max_gaplength_to_model 8 \
	-MR::fast \
	-nstruct 1 \
	-ignore_unrecognized_res -overwrite
