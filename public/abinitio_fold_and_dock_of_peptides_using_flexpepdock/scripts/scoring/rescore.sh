#!/bin/csh

#given a score file that is the result of a fpdock-abinitio run, adds as the last column 
#the 'reweighted score' 

set scorefile=$1

set sc=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="score") print i}}'`
set I_sc=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="I_sc") print i}}'`
set pep_sc=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="pep_sc") print i}}'`

cat $scorefile | awk '{if (NR==1) {print $0,"reweighted_score"} else {print $0,$'$sc'+$'$I_sc'+'$pep_sc'}}' > newscore.sc

echo Rescoring done - created newscore.sc
echo
