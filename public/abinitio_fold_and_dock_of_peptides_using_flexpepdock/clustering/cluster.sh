#!/bin/bash

printScoreFile_byHeader.pl ../../score.sc I_sc reweighted_sc description | sort -nk 3 | head -500 | awk '{print $NF}' >pdb_list

PATH_TO_EXE="rosetta/main/source/bin"
PATH_TO_DB="rosetta/main/database"

radius=$1
native=$2
decoys=$3

#calculate actual clustering radius taking intoaccount protein lenght.
len=`grep CA $native | wc -l`
plen=`awk '$5=="B"' $native | grep CA | wc -l`
actualR=`date | awk '{print sqrt('$plen'/'$len')*'$radius'}'`
echo actual radius is $actualR

#send clustering run
if [ -e c.0.0.pdb ]; then
	rm c.*.pdb
fi

if [ ! -e cluster.silent ]; then
		$PATH_TO_EXE/cluster.linuxgccrelease -in:file:silent $decoys -in:file:silent_struct_type binary -database $PATH_TO_DB -cluster:radius $actualR -in:file:fullatom -tags `cat pdb_list` -silent_read_through_errors > clog
    echo Done clustering.
fi

echo Now printing results.

x=`wc -l pdb_list | awk '{print $1}'`
tail -$((x+5)) clog | head -${x} | awk '{print $4,$5,$6}' > cluster_list

if [ -e pdb_list_sc ];then rm pdb_list_sc;fi
for i in `awk '{print $1}' cluster_list`;do sed 1d ../../score.sc | head -1 > tmp1; grep $i ../../score.sc >> tmp1; printScoreFile_byHeader.pl tmp1 I_sc pep_sc reweighted_sc rmsBB_if bestRMS_4mer_all bestRMS_6mer_all rmsALL_if | tail -1 | awk '{print $2,$3,$4,$5,$6,$7,$8}' >>pdb_list_sc; done
paste cluster_list pdb_list_sc >cluster_list_sc

sort -nk 4 cluster_list_sc | sort -u -k2,2 | sort -nk 4 | head -10 >cluster_list_I_sc_sorted 
sort -nk 5 cluster_list_sc | sort -u -k2,2 | sort -nk 5 | head -10 >cluster_list_pep_sc_sorted 
sort -nk 6 cluster_list_sc | sort -u -k2,2 | sort -nk 6 | head -10 >cluster_list_reweighted_sc_sorted 


exit 0
