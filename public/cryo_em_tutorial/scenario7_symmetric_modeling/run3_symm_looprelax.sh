#!/bin/bash

$ROSETTA3_SRC/bin/loopmodel.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in::file::centroid_input \
	-symmetry:symmetry_definition 1K4C_edited.symm \
	-symmetry::initialize_rigid_body_dofs \
	-in::file::s 1K4C_mem_abrelax_INPUT.pdb \
	-loops:loop_file 1K4C.loopfile \
	-loops::frag_sizes 9 3 1 \
	-loops::frag_files aa1kcs_09_05.200_v1_3  aa1kcs_03_05.200_v1_3 none \
	-loops::remodel quick_ccd \
	-loops::intermedrelax no \
	-loops::refine no \
	-loops::relax relax \
	-relax::default_repeats 1 \
	-relax::jump_move true \
	-score::weights score12_nosol.wts \
	-loops::build_attempts 10 \
	-edensity::mapfile 1K4C.5A.mrc \
	-edensity::mapreso 8.0 \
	-edensity::grid_spacing 4.0 \
	-edensity::whole_structure_ca_wt 0.05 \
	-edensity::score_symm_complex true \
	-overwrite
