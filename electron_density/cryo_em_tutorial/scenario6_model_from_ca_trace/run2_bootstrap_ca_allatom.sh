#!/bin/bash

$ROSETTA3_SRC/bin/ca_to_allatom.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in:file:s 1bbh_threaded.pdb \
	-in:file:fullatom \
	-in:file:native 1bbh.pdb \
	-edensity::mapfile 1bbh.5A.mrc \
	-edensity::realign min \
	-edensity::mapreso 5.0 \
	-edensity::grid_spacing 2.0 \
	-edensity::atom_mask 3.5 \
	-edensity::whole_structure_allatom_wt 0.05 \
	-RBSegmentRelax::cst_wt 0.1 \
	-RBSegmentRelax::cst_width 2.0 \
	-RBSegmentRelax::rb_scorefxn score12_full \
	-RBSegmentRelax::rb_file 1bbh.rbsegs \
	-RBSegmentRelax::nrbmoves 20 \
	-RBSegmentRelax::helical_movement_params 30.0 0.5 2.0 0.5 \
	-loops::vall_file /scratch/ROSETTA/nnmake_database/vall.dat.2006-05-05 \
	-loops::frag_sizes 9 3 1 \
	-loops::frag_files aa1bbh_09_05.200_v1_3.gz aa1bbh_03_05.200_v1_3.gz none \
	-loops::random_loop \
	-loops::build_attempts 10 \
	-loops::remodel quick_ccd \
	-loops::intermedrelax no \
	-loops::refine no \
	-loops::relax relax \
	-relax::default_repeats 1 \
	-out::suffix _allatom \
	-overwrite
