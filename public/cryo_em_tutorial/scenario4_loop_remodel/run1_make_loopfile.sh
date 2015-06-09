#!/bin/bash

$ROSETTA3_SRC/bin/loops_from_density.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in::file::s 1cid_threaded.pdb \
	-in:file:fullatom \
	-edensity::mapfile 1cid.5A.mrc \
	-edensity::whole_structure_allatom_wt 1.0 \
	-edensity::realign min \
	-edensity::sliding_window 9 \
	-edensity::mapreso 4 \
	-edensity::grid_spacing 2 \
	-max_helix_melt -1 \
	-max_strand_melt 3 \
	-frac_loop 0.2

