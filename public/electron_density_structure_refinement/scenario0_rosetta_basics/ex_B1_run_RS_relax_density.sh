#!/bin/bash

~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::s 1isrA.pdb \
 -parser::protocol ex_B1_run_RS_relax_density.xml \
 -ignore_unrecognized_res \
 -edensity::mapreso 5.0 \
 -edensity::cryoem_scatterers \
 -crystal_refine \
 -out::suffix _relax \
 -default_max_cycles 200

