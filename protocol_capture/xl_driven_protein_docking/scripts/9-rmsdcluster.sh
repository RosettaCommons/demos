#!/bin/bash
if [ $# -ne 4 ]
then
    echo
    echo "Usage: $0 path2cleftXplorer.jar path2colt.jar path2cdk-1.0.2.jar path2R"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    localClusterDir="$projectDir/cluster/local"
    bestModelsDir="$projectDir/bestModels";
    rscript="$pwd/scripts/rmsd_cluster.R"

    allVsAllFile="$localClusterDir/rmsd_allVsAll.txt"
    clusterFile="$localClusterDir/cluster.txt"
    bestModelsFile="$localClusterDir/bestModels.txt"

    echo
    cd $localClusterDir
    
    # run Superimpose
    receptorChainId=`grep -h ^ATOM *.pdb | head -1 | perl -ne '$c=substr($_, 21, 1); print $c'`
    cmd="java -Xmx512m -cp $1:$2:$3 cX/Superimpose -dir ./ -chain $receptorChainId:$receptorChainId -xseq -ca -r -dock"
    echo "Calculating all vs all ligand RMSD values: $cmd"
    eval $cmd > $allVsAllFile

    cmd="$4 CMD BATCH $rscript"
    echo "Computing a hierarchical clustering with complete linkage: $cmd"
    eval $cmd


    echo "Selecting following lowest scoring models from largest 3 clusters as best predictions."
    perl -ane '$s=`grep ^pose $F[0] | sed "s/.* //"`; $n{$F[1]}++; if(exists $c{$F[1]}){if($s<$c{$F[1]}){$c{$F[1]}=$s;$d{$F[1]}=$F[0]}}else{$c{$F[1]}=$s;$d{$F[1]}=$F[0]; END{foreach my $c (sort {$b<=>$a} keys %n){print "$d{$c}\t$c{$c}\t$n{$c}\t$c\n"}}}' $clusterFile > $bestModelsFile

    echo "Copying following top3 cluster representative to $bestModelsDir:"
    mkdir $bestModelsDir
    for p in `head -3 $bestModelsFile | cut -f1`; do 
	echo $p
	cp $p $bestModelsDir
    done

    echo "DONE!"
    echo
    cd $pwd
fi
