#!/bin/bash
if [ $# -ne 2 ]
then
    echo
    echo "Usage: $0 path2XwalkBinDir path2naccess"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    localClusterDir="$projectDir/cluster/local"
    
    clusterFilesArchive="allClusterFiles.tgz"
    bestModelDir="$projectDir/bestModels"
    interfaceFile="$bestModelDir/predicted_interface.pdb"

    echo
    cd $localClusterDir
    
    # create an archive of all PDB files having remained in the cluster directory, indicating that these structures satisfy to most cross-links and have a sufficiently large binding interface.
    tar cfz $clusterFilesArchive *.pdb

    # run AvgInterface
    cmd="java -cp $1 AvgInterface $clusterFilesArchive `ls $bestModelDir/*.pdb | head -1` $2"
    echo "Computing contact frequencies with AvgInterface and saving resulting PDB file in $interfaceFile: $cmd"
    eval $cmd > $interfaceFile

    echo "DONE!"
    echo
    cd $pwd
fi
