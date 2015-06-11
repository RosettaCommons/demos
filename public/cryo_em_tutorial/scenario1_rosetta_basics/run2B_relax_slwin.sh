#!/bin/bash

$ROSETTA3_SRC/bin/relax.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s 2JEL_P_idl.pdb \
 -relax::fast \
 -relax::default_repeats 4 \
 -relax::jump_move true \
 -edensity::mapfile 1POH.5A.mrc \
 -edensity::mapreso 5.0 \
 -edensity::grid_spacing 2.0 \
 -edensity::sliding_window_wt 0.5 \
 -edensity::sliding_window 3 \
 -out:suffix _slwin \
 -overwrite
