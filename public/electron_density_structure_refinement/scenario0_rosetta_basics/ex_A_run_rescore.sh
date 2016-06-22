#!/bin/bash

~/Rosetta/main/source/bin/score_jd2.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::s 1isrA.pdb 1issA.pdb \
 -ignore_unrecognized_res \
 -edensity::mapfile 1issA_6A.mrc \
 -edensity::mapreso 5.0 \
 -edensity::grid_spacing 2.0 \
 -edensity::fastdens_wt 20.0 \
 -edensity::cryoem_scatterers \
 -crystal_refine
