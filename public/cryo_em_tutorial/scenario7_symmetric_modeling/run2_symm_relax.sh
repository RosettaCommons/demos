#!/bin/bash

$ROSETTA3_SRC/bin/relax.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in:file:s 1K4C_mem_abrelax_INPUT.pdb \
	-symmetry:symmetry_definition 1K4C.symm \
	-symmetry::initialize_rigid_body_dofs \
	-score::weights score12_nosol.wts \
	-relax::default_repeats 1 \
	-relax::jump_move true \
	-edensity::mapfile 1K4C.5A.mrc \
	-edensity::mapreso 5.0 \
	-edensity::grid_spacing 2.0 \
	-edensity::whole_structure_allatom_wt 0.1 \
	-edensity::score_symm_complex true \
	-overwrite
