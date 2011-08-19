#!/bin/bash


# This script converts restraints in XPlor format into Rosetta's cst format.
# ( block_276112.mr -> gb1_tedor.cst in this particular case)
# To run it, download bioshell distribution from:
#	http://bioshell.chem.uw.edu.pl
# and put the jars on your CLASSPATH

XPLOR_FILE=$1
PDB_FILE=$2

java apps.Restraints -distances_xplor=$XPLOR_FILE -ip=$PDB_FILE  -select.bb -output_restraints=all.rst
tr -d '_' < all.rst | awk '{if($3!=$6) print "AtomPair",$4,$3,$7,$6,"LINEAR_PENALTY", $10,0.0,2/$11,1.0}' > out.cst
rm all.rst