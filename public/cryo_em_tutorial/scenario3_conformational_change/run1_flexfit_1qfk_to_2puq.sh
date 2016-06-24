#!/bin/bash

~/Rosetta/main/source/bin/relax.default.linuxclangrelease \
 -s 1QFK_L.pdb \
 -relax::default_repeats 1 \
 -relax::jump_move true \
 -edensity::mapfile 2PUQ.5A.mrc \
 -edensity::mapreso 12.0 \
 -edensity::grid_spacing 6.0 \
 -edensity::ca_mask 20.0 \
 -edensity::whole_structure_ca_wt 0.1 \
 -overwrite

mv 1QFK_L_0001.pdb 1QFK_L_lowres_fit.pdb
