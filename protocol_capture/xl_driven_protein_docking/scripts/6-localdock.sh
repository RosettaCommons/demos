#!/bin/bash
if [ $# -ne 5 ]
then
    echo
    echo "Usage: $0 path2docking_protocol path2rosettaDatabase numberOfModels nativeFile cstFile"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    localDockingDir="$projectDir/docking/local"

    flagFile="$pwd/inputs/local_docking.flags"
    modelsFile="$localDockingDir/docking-models.tgz"    

    echo
    cd $localDockingDir

    chainIds=`grep -h ^ATOM *.pdb | perl -ne '$c=substr($_,21,1); if($chains!~/$c/){$chains.=$c."_"} END{chop($chains); print "$chains"}'`

    # do local docking refinement
    for p in *.pdb; do
	cmd="$1 -database $2 -in:file:s $p -docking:partners $chainIds -in:file:native $pwd/$4 -constraints:cst_file $pwd/$5 -out:nstruct $3 @$flagFile"
	echo "Running local refinement docking calculation: $cmd"
	eval $cmd
    done

    echo "Storing PDB files in an archive before deleting them. This will relieve the file system"
    tar cfz $modelsFile *_????_????.pdb && rm -f *_????_????.pdb
    
    echo "DONE!"
    echo
    cd $pwd
fi
