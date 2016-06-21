#!/bin/bash

XPLOR_FILE=$1
PDB_FILE=$2

java apps.Restraints -distances_xplor=$XPLOR_FILE -ip=$PDB_FILE  -select.bb -output_restraints=all.rst
tr -d '_' < all.rst | awk '{if($3!=$6) print "AtomPair",$4,$3,$7,$6,"LINEAR_PENALTY", $10,0.0,2/$11,1.0}' > out.cst
rm all.rst