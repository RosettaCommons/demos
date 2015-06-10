#!/bin/bash
if [ $# -ne 4 ]
then
    echo
    echo "Usage: $0 path2XwalkBinDir dockingDir clusterDir2beCreated xwalkXLsFile"
    echo
else
    SIFS=$IFS;
    IFS=$'\n';

    pwd=`pwd`
    projectDir="$pwd/output"

    modelsFile="$pwd/$2/docking-models.tgz"
    distsFile="$pwd/$2/docking-models_dist.txt"
    satisfyFile="$pwd/$2/conform.txt"

    echo
    cd $2
    
    # run Xwalk
    cmd="java -Xmx512m -cp $1 Xwalk -infile $modelsFile -dist $pwd/$4 -max 34 -xSC -mono -out $distsFile -f"
    echo "Running Xwalk: $cmd"
    eval $cmd

    echo "Determining PDB files that satisfy to most cross- and mono-links." 
    perl -ne '@a=split(/\t/); if(($a[6]>0 and $a[6]<=34) or ($#a==3 and $a[3]==1)){$h{$a[1]}++; if($h{$a[1]}>$max){$max=$h{$a[1]}}}; END{foreach $n(keys %h){print "$n\t$h{$n}\n" if($h{$n}==$max)}}' $distsFile > $satisfyFile


    echo "Select the top 500 structures that all satisfy to `cut -f2 $satisfyFile | head -1` cross- and mono-links and have the lowest ROSETTA energy score:"
    cd $pwd
    mkdir -p $3
    cd $3

    i=0
    for s in `grep SCORE $pwd/$2/score*sc | grep -v description | sort -k2 -n `; do
	i=$(($i+1))
	m=`echo $s | sed "s/.* //"`;
	if [ `grep -c $m.pdb $satisfyFile` -ne 0 -a $i -le 500 ] 
	then
	    echo "$i. $m.pdb"
	    tar xfz $modelsFile $m.pdb
	fi
    done

    echo "DONE!"
    echo
    cd $pwd

    IFS=$SIFS;
fi
