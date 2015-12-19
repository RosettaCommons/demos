#!/bin/bash

rm default.sc

$ROSETTA3_SRC/bin/score.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s 2JEL_P.pdb 2JEL_P_idl.pdb 2JEL_P_idl_0001.pdb 2JEL_P_idl_cst_0001.pdb 2JEL_P_idl_slwin_0001.pdb 1POH_A.pdb \
 -in::file::native 1POH_A.pdb \
 -score:weights score13_env_hb \
 -ignore_unrecognized_res \
 -edensity::mapfile 1POH.5A.mrc \
 -edensity::mapreso 5.0 \
 -edensity::grid_spacing 2.0 \
 -edensity::whole_structure_allatom_wt 0.1
