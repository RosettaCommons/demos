#!/bin/csh

# assumes gnuin is present in running dir

# usage: cluster.sh <pdb name> <top x structures> <radius> <score.sc> <native.pdb> <decoys.silent> <stype=column number of score according to which to select the top x structures>
# e.g. cluster.sh 2FMF 100 1 ../score.sc ../native.pdb ../decoys.silent 2 ---> clusters top 100 by total score with 1A peptide bb rms radius.

set PATH_TO_EXE=/vol/ek/ravehb/rosetta/svn_mini/mini/bin
set PATH_TO_DB=/vol/ek/ravehb/rosetta/svn_mini/minirosetta_database/

set name=$1
set top=$2
set radius=$3
set score=$4
set native=$5
set decoys=$6
set stype=$7


#calculate actual clustering radius taking intoaccount protein lenght.
set len=`grep CA $native | wc -l`
set plen=`awk '$5=="B"' $native | grep CA | wc -l`
set actualR=`date | awk '{print sqrt('$plen'/'$len')*'$radius'}'`
#send clustering run
if (! -e c.0.0.pdb ) then
	$PATH_TO_EXE/cluster.linuxgccrelease -in:file:silent $decoys -in:file:silent_struct_type binary -database $PATH_TO_DB -cluster:radius $actualR -in:file:fullatom -tags `cat $score | sort -nrk $stype | tail -$top | awk '{print $(NF-1)}'` -silent_read_through_errors >! clog
endif
#find clusters in scorefile
#cat clog | awk 'BEGIN{s=0}{if ($0 ~ /Sorting/) s=1; if (s==1) print $0}' | awk 'BEGIN{c=-1}{if ($3 ~ /Cluster:/) c++; else {if (NR>2) print "echo -n `grep "$3" '$score'`; echo \" "c"\"" } if (c==10) exit;}' > ! getFromSc.tmp

#cat clog | awk 'BEGIN{s=0}{if ($0 ~ /Summary/) s=1; if (s==1) print $0}' | tail -n +2 | awk 'BEGIN {end=0; count=0; sum=0;} {if ( $3 ~ /------/) end=1; if ( end==0 ) {print $0, $4} }'

set topW5=0
@ topW5 = $top + 5

#join by score
tail -$topW5 clog | head -$top | awk '{print $4, $5 }' | sort -nk2 | sort >! tmp.listA
cat newscore.sc | awk '{print $(NF-1),$NF}' | sort > tmp.listB
echo description cluster reweight_sc
join -1 1 -2 1 tmp.listA tmp.listB | sort -nk2 -k 3 | awk 'BEGIN{cur=0}{if(cur==$2){print; cur++;}}' | sort -nk 3

rm -rf tmp.list?
