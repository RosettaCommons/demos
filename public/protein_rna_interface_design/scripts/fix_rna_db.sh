#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 /path/to/minirosetta_database"
	exit
fi

db=$1

cp ${db}/chemical/residue_type_sets/rna/residue_types/*.params ${db}/chemical/residue_type_sets/fa_standard/residue_types/nucleic

cp ${db}/chemical/residue_type_sets/rna/patches/LowerRNA.txt ${db}/chemical/residue_type_sets/fa_standard/patches
cp ${db}/chemical/residue_type_sets/rna/patches/UpperRNA.txt ${db}/chemical/residue_type_sets/fa_standard/patches

perl -ni -e 'if ($_ =~ /\#\# Nucleic Acid Types/) { print $_; print "residue_types/nucleic/RAD.params\nresidue_types/nucleic/RCY.params\nresidue_types/nucleic/RGU.params\nresidue_types/nucleic/URA.params\n"} else { print $_ }' ${db}/chemical/residue_type_sets/fa_standard/residue_types.txt
sed -i '$a\\n## RNA patch files\npatches/LowerRNA.txt\npatches/UpperRNA.txt' ${db}/chemical/residue_type_sets/fa_standard/patches.txt