#!/bin/bash

$ROSETTA3_SRC/bin/loopmodel.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in::file::s 1cid_threaded.pdb \
	-loops::loop_file 1cid_threaded.loopfile \
	-nstruct 1 \
	-in::file::fullatom \
	-loops::remodel quick_ccd_moves \
	-loops::intermedrelax no \
	-loops::refine no \
	-loops::relax relax \
	-loops::frag_sizes 9 3 1 \
	-loops::frag_files aa1cid_09_05.200_v1_3.gz  aa1cid_03_05.200_v1_3.gz none \
	-loops::extended \
	-loops::random_grow_loops_by 4 \
	-relax::default_repeats 1 \
	-relax::jump_move true \
	-edensity::mapfile 1cid.5A.mrc \
	-edensity::realign min \
	-edensity::mapreso 5.0 \
	-edensity::grid_spacing 2.5 \
	-edensity::whole_structure_allatom_wt 0.05 \
	-overwrite
