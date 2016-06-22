#!/bin/bash

# get chain A of input only
outfile=`echo $1 | sed 's/\.pdb/_chainA.pdb/'`
cat $1 | awk '{ if ( "A"==substr($0,22,1)) print $0 }' > $outfile

# run the script
~/rosetta_workshop/rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::s $outfile \
 -parser::protocol ex_C1_bfact_FSC.xml \
 -ignore_unrecognized_res \
 -edensity::mapreso 3.4 \
 -edensity::cryoem_scatterers \
 -out::suffix _bfact \
 -crystal_refine
