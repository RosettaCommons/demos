#!/bin/bash

~/rosetta_workshop/rosetta/main/source/bin/density_tools.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::s 3j5p_transmem.pdb \
 -edensity::mapfile TRPV1_half2.mrc \
 -edensity::mapreso 3.4 \
 -edensity::cryoem_scatterers \
 -crystal_refine \
 -denstools::verbose

