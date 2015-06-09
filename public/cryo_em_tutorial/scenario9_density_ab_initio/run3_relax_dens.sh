#!/bin/bash


$ROSETTA3_SRC/bin/relax.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s S_0001_idl.pdb \
 -relax::default_repeats 1 \
 -relax::jump_move true \
 -edensity::mapfile 1nsf.5A.mrc \
 -edensity::mapreso 5.0 \
 -edensity::grid_spacing 2.0 \
 -edensity::whole_structure_ca_wt 0.1 \
 -overwrite
