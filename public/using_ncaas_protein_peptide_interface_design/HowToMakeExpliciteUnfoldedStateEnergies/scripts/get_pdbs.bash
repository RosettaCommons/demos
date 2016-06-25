#!/bin/bash

# get_pdbs.bash - a simple script to get pdbs from a PICIES culled pdb formated file

# get list of four letter codes form file

curl -O http://dunbrack.fccc.edu/Guoli/culledpdb_hh/$1.gz
gunzip $1

flcs=`tail +2 $1 | cut -c 1-4 | tr '[:upper:]' '[:lower:]' | while read line ; do curl -O http://files.rcsb.org/view/$line.pdb  ; done`


 
