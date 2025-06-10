#!/bin/bash
# cascaded_prepare.sh <ligand name> <number of top ligands> <protein name>
# This script prepares the inputs for cascaded sampling. Run this script under the production folder

currLIG=$1
nTOP=$2
protein=$3

cd $currLIG/
Top=($( tail -n +2 design_scores.sc | head -n $nTOP | awk '{print $1}' ))
INDEX=1
now=$(date)
echo "=====$now=====" >> ../Cascaded_README
echo "Current ligand: $currLIG" >> ../Cascaded_README
echo "Selected $nTOP ligands for cascaded sampling." >> ../Cascaded_README

for lig in "${Top[@]}"
do
	complex="$lig.pdb"
	newLIG="${lig}_ligand.pdb"
	grep -- "^HETATM" "$complex" > "$newLIG"
	cp $newLIG ../../inputs/${currLIG}_LIG$INDEX.pdb
	title=$( awk 'NR==1{print $4}' $newLIG )
	obabel ../../inputs/${currLIG}_LIG$INDEX.pdb -O ../../inputs/${currLIG}_LIG$INDEX.sdf --title $title
	obabel ../../inputs/${currLIG}_LIG$INDEX.sdf -O ../../inputs/${currLIG}_LIG$INDEX.png --gen3D -d --title
	cat ../../inputs/$protein.pdb ../../inputs/${currLIG}_LIG$INDEX.pdb > ../../inputs/${protein}_${currLIG}_LIG$INDEX.pdb
	echo "$lig -> ${currLIG}_LIG$INDEX" >> ../Cascaded_README
	((INDEX = INDEX +1))
done
