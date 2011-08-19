#!/bin/bash

cd generation2

$ROSETTA3_SRC/bin/loops_from_density.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in::file::s *.pdb \
	-in:file:fullatom \
	-edensity::mapfile ../1q0p.5A.mrc \
	-edensity::whole_structure_allatom_wt 1.0 \
	-edensity::realign min \
	-edensity::sliding_window 9 \
	-edensity::mapreso 8 \
	-edensity::grid_spacing 4 \
	-max_helix_melt -1 \
	-max_strand_melt 2 \
	-frac_loop 0.4

cd ..
