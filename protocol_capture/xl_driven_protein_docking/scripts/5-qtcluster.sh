#!/bin/bash
if [ $# -ne 3 ]
then
    echo
    echo "Usage: $0 path2cleftXplorer.jar path2colt.jar path2cdk-1.0.2.jar"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    globalClusterDir="$projectDir/cluster/global"
    localDockingDir="$projectDir/docking/local"
    qtPl="$pwd/scripts/qt_cluster.pl"

    allVsAllFile="$globalClusterDir/qt_allVsAll.txt"
    clusterFile="$globalClusterDir/cluster.txt"

    echo
    cd $globalClusterDir
    
    # run Superimpose
    receptorChainId=`grep -h ^ATOM *.pdb | head -1 | perl -ne '$c=substr($_, 21, 1); print $c'`
    cmd="java -Xmx512m -cp $1:$2:$3 cX/Superimpose -dir ./ -chain $receptorChainId:$receptorChainId -dock -xseq -ca -ta 3:8"
    echo "Calculating all vs all quality thresholds (QT): $cmd"
    eval $cmd > $allVsAllFile

    cmd="perl $qtPl -in $allVsAllFile -dir ./ -v"
    echo "Determining three largest QT clusters: $cmd"
    eval $cmd > $clusterFile

    echo "copying following cluster representatives to $localDockingDir:"
    mkdir -p $localDockingDir
    for p in `head -3 $clusterFile | cut -f1`; do 
	echo $p
	cp $p $localDockingDir
    done

    echo "DONE!"
    echo
    cd $pwd
fi
