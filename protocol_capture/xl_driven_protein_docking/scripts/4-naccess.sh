#!/bin/bash
if [ $# -ne 3 ]
then
    echo
    echo "Usage: $0 path2naccess removeFilesWith2smallInterface clusterDir"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    bsaPl="$pwd/scripts/bsa.pl"

    echo
    cd $3
    
    echo "Calculating interface sizes with NACCESS"
    for p in *.pdb; do
	bsa=`perl $bsaPl $p -naccess $1 | cut -f 4`
	echo -n "$p has a binding surface area of $bsa"
	perl -e 'if('$bsa' < 900){print ", which is smaller than 900. "; if('$2'){print "Removing it!"; unlink("'$p'")}else{print "Should be removed!"}}';
	echo
    done 
    
    echo "DONE!"
    echo
    cd $pwd
fi
