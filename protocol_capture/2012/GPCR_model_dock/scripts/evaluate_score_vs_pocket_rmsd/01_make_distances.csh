#!/bin/tcsh

set wdr = `pwd`

foreach i (`ls -1 structures/`) # ls -1 structures/
	ls -1 structures/$i/ | sed 's/.pdb//g' > pdb_lists/$i.lst
	sed 's/__pdbidChain__/'$i'/g' scripts/01_distances.py |sed 's,__working_dir__,'$wdr/',g' > scripts/01_distances_$i.py
	pymol -c scripts/01_distances_$i.py > & scripts/01_distances_$i.log
end
