#!/bin/bash

$ROSETTA3_SRC/bin/mr_protocols.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-MR::mode cm \
	-in::file::extended_pose 1 \
	-in::file::fasta 1XVQ.fasta \
	-in::file::alignment 1xvq_2bmx.fasta \
	-in::file::template_pdb 2bmxA.pdb \
	-cm::min_loop_size 4 \
	-cm::loop_close_level 0 \
	-cm::loop_rebuild_filter 150 \
	-MR::max_gaplength_to_model 999 \
	-loops::frag_sizes 9 3 1 \
	-loops::frag_files aaxvqn_09_05.200_v1_3.gz aaxvqn_03_05.200_v1_3.gz none \
	-loops::remodel quick_ccd \
	-loops::relax relax \
	-relax::default_repeats 1 \
	-relax::jump_move true \
	-edensity::mapfile 1XVQ.5A.mrc \
	-edensity::realign min \
	-edensity::mapreso 5.0 \
	-edensity::grid_spacing 2.5 \
	-edensity::whole_structure_allatom_wt 0.1 \
	-overwrite

