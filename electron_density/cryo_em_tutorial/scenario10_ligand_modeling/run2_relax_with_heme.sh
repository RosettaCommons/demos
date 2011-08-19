#!/bin/bash

$ROSETTA3_SRC/bin/relax.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s 1bbh_threaded_with_heme.pdb \
 -relax::default_repeats 1 \
 -relax::jump_move true \
 -extra_res_fa HEM.params \
 -edensity::mapfile 1bbh.5A.mrc \
 -edensity::mapreso 5.0 \
 -edensity::grid_spacing 2.5 \
 -edensity::whole_structure_allatom_wt 0.1
