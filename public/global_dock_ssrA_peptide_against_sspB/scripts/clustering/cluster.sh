#!/bin/csh

# assumes gnuin is present in running dir

# usage: cluster.sh <top x structures> <radius> <score.sc> <native.pdb> <decoys.silent> <stype=column number of score according to which to select the top x structures>
# e.g. cluster.sh 100 1 ../score.sc ../native.pdb ../decoys.silent 2 ---> clusters top 100 by total score with 1A peptide bb rms radius.

set PATH_TO_EXE=/vol/ek/ravehb/rosetta/svn_mini/rosetta_source/bin
set PATH_TO_DB=/vol/ek/ravehb/rosetta/svn_mini/rosetta_database/
set PATH_TO_SCRIPTS=/vol/ek/ravehb/rosetta/RosettaCon2011/demos/global_dock_ssrA_peptide_against_sspB/scripts/
set pepChain=B

set top=$1
set radius=$2
set scorefile=$3
set native=$4
set decoys_file=$5
set scoretype=$6


#calculate actual clustering radius taking into account protein lenght.
echo BEGIN CLUSTERING of $scorefile, decoys file $decoys_file, using score-type $scoretype
echo radius $radius top $top, native $native
set len=`awk '/^ATOM/ && substr($0,14,3)=="CA "' $native | wc -l`
set plen=`awk '/^ATOM/ && substr($0,14,3)=="CA " && substr($0,22,1)=="'"$pepChain"'"' $native | wc -l`
set actualR=`date | awk '{print sqrt('$plen'/'$len')*'$radius'}'`
echo Peptide clustering radius $radius A 
echo "(actual radius adjusted to $actualR A for total-length=$len and peptide-length=$plen)"
#send clustering run
if (! -e c.0.0.pdb ) then
	$PATH_TO_EXE/cluster.linuxgccrelease -in:file:silent $decoys_file -in:file:silent_struct_type binary -database $PATH_TO_DB -cluster:radius $actualR -in:file:fullatom -tags `$PATH_TO_SCRIPTS/utils/printScoreFile_byHeader.pl $scorefile $scoretype description | sort -nrk 2 | tail -$top | awk '{print $3}'` -silent_read_through_errors >! clog
endif
#find clusters in scorefile
#cat clog | awk 'BEGIN{s=0}{if ($0 ~ /Sorting/) s=1; if (s==1) print $0}' | awk 'BEGIN{c=-1}{if ($3 ~ /Cluster:/) c++; else {if (NR>2) print "echo -n `grep "$3" '$score'`; echo \" "c"\"" } if (c==10) exit;}' > ! getFromSc.tmp

#cat clog | awk 'BEGIN{s=0}{if ($0 ~ /Summary/) s=1; if (s==1) print $0}' | tail -n +2 | awk 'BEGIN {end=0; count=0; sum=0;} {if ( $3 ~ /------/) end=1; if ( end==0 ) {print $0, $4} }'

# report top clusters ranked by $scoretype
set topW5=0
@ topW5 = $top + 5
tail -$topW5 clog | head -$top | awk '{print $4, $5 }' | sort -nk2 | sort >! tmp.listA
$PATH_TO_SCRIPTS/utils/printScoreFile_byHeader.pl $scorefile description $scoretype | sort > tmp.listB
echo description cluster $scoretype
join -1 1 -2 1 tmp.listA tmp.listB | sort -nk2 -k 3 | awk 'BEGIN{cur=0}{if(cur==$2){print; cur++;}}' | sort -nk 3 > cluster_results.txt
rm -rf tmp.list?
