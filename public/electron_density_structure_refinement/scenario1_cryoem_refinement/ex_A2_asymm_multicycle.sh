#!/bin/bash

~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::s 3j5p_transmem_A.pdb \
 -parser::protocol ex_A2_asymm_multicycle.xml \
 -ignore_unrecognized_res \
 -nstruct 5 \
 -edensity::mapreso 3.4 \
 -edensity::cryoem_scatterers \
 -in:file:centroid_input \
 -out::suffix _asymm_multicyc \
 -crystal_refine
