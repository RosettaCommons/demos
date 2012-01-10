#!/bin/bash

$ROSETTA3_SRC/bin/relax.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s 2PUQ_L.pdb \
 -relax::default_repeats 1 \
 -relax::jump_move true \
 -edensity::mapfile 1QFK.5A.mrc \
 -edensity::mapreso 12.0 \
 -edensity::grid_spacing 6.0 \
 -edensity::ca_mask 20.0 \
 -edensity::whole_structure_ca_wt 0.2 \
 -overwrite

mv 2PUQ_L_0001.pdb 2PUQ_L_lowres_fit.pdb