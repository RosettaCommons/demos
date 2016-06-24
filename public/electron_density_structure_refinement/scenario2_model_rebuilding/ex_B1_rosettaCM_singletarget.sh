#!/bin/bash

~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in:file:fasta t20s.fasta \
 -parser:protocol ex_B1_rosettaCM_singletarget.xml \
 -nstruct 50 \
 -relax:minimize_bond_angles \
 -relax:min_type lbfgs_armijo_nonmonotone \
 -relax:jump_move true \
 -relax:default_repeats 2 \
 -out::suffix _singletgt \
 -edensity::mapfile t20S_41A_half1.mrc \
 -edensity::mapreso 5.0 \
 -edensity::cryoem_scatterers \
 -default_max_cycles 200
