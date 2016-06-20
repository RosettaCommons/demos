#!/bin/bash
if [ $# -ne 5 ]
then
    echo
    echo "Usage: $0 path2docking_protocol path2rosettaDatabase numberOfModels nativeFile cstFile"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    globalDockingDir="$projectDir/docking/global"

    inputFile="$projectDir/relax/complex-relaxed.pdb"
    flagFile="$pwd/inputs/global_docking.flags"
    modelsFile="$globalDockingDir/docking-models.tgz"    

    echo
    mkdir -p $globalDockingDir
    cd $globalDockingDir

    chainIds=`grep -h ^ATOM $inputFile | perl -ne '$c=substr($_,21,1); if($chains!~/$c/){$chains.=$c."_"} END{chop($chains); print "$chains"}'`

    # do global low resolution docking
    cmd="$1 -database $2 -in:file:s $inputFile -docking:partners $chainIds -in:file:native $pwd/$4 -constraints:cst_file $pwd/$5 -out:nstruct $3 @$flagFile"
    echo "Running global docking calculation: $cmd"
    eval $cmd

    echo "Storing PDB files in an archive before deleting them. This will relieve the file system"
    tar cfz $modelsFile *.pdb && rm -f *.pdb
    
    echo "DONE!"
    echo
    cd $pwd
fi
