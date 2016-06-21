#!/bin/bash

~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::s 3j5p_transmem_A.pdb \
 -parser::protocol ex_A1_asymm_1cycle_legacy.xml \
 -ignore_unrecognized_res \
 -nstruct 1 \
 -edensity::mapreso 3.4 \
 -edensity::cryoem_scatterers \
 -in:file:centroid_input \
 -out::suffix _refine_1cyc \
 -default_max_cycles 200 \
 -crystal_refine
