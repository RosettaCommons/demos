#!/bin/bash

# get_pdbs.bash - a simple script to get pdbs from a PICIES culled pdb formated file

# ftp host for the pdbs
HOST=ftp://ftp.wwpdb.org

# get list of four letter codes form file
flcs=`tail +2 $1 | cut -c 1-4 | tr '[:upper:]' '[:lower:]' | while read line ; do echo "get pdb$line.ent.gz" ; done`

# redirect it all in ftp to download pbds
ftp -iv $HOST << EOF
cd pub/pdb/data/structures/all/pdb/
$flcs 
bye
EOF  
