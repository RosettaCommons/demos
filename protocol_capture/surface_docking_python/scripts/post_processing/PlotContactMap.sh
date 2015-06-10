#!/bin/bash

echo "Generating contact map"
ls *.pdb > PDBsInDir
DecoyCnt=$( wc -l PDBsInDir | cut -f1 -d" " )

Current_structure=1
while [ $Current_structure -le $DecoyCnt ]
  do
  Current_structure_name=$( head -$Current_structure PDBsInDir | tail -1 )
  grep ATOM $Current_structure_name > protein.pdb
  contacts.py protein.pdb > $Current_structure_name.con
  rm protein.pdb
  Current_structure=$[ $Current_structure + 1 ] 
done

cat AdsState*.con > Ads.con
cat SolState*.con > Sol.con

# Should work if there's either one or two types of files
ads=$( find -empty -name Ads.con )
#sol=$( find -empty -name Sol.con )

if [ -n "$ads" ]
then 
    ContactMap.py Sol.con native.pdb.con 
else
    ContactMap.py Ads.con native.pdb.con 
fi

rm *con PDBsInDir
    
~/Applications/bin/gnuplot ContactMap.gnuplot

echo "Contact map DONE"
