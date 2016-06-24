#!/bin/bash
if [ $# -ne 5 ]
then
    echo
    echo "Usage: $0 receptor.pdb ligand.pdb path2relaxOrprepack path2rosettaDatabase doRelax"
    echo
else
    pwd=`pwd`
    projectDir="$pwd/output"
    relaxDir="$projectDir/relax"

    complexFile="$relaxDir/complex.pdb"
    relaxedComplex="$relaxDir/complex-relaxed.pdb"
    relaxedReceptorFile="$relaxDir/"`echo $1 | sed "s/.*\///" | sed "s/.pdb/_0001.pdb/"`
    relaxedLigandFile="$relaxDir/"`echo $2 | sed "s/.*\///" | sed "s/.pdb/_0001.pdb/"`

    echo
    mkdir -p $relaxDir
    cd $relaxDir

    receptorChainId=`grep -h ^ATOM $pwd/$1 | head -1 | perl -ne '$c=substr($_, 21, 1); $c="A" if($c eq " "); print $c'`
    ligandChainId=`grep -h ^ATOM $pwd/$2 | head -1 | perl -ne '$c=substr($_, 21, 1); $c="B" if($c eq " "); print $c'`
    if [ $5 -gt 0 ]; 
	then
        # do relax
	cmd="$3 -in:file:s $pwd/$1 -out:nstruct 1 -database $4 -relax:sequence -constrain_relax_to_start_coords"
	echo "Relaxing receptor structure: $cmd"
	eval $cmd
	cmd="$3 -in:file:s $pwd/$2 -out:nstruct 1 -database $4 -relax:sequence -constrain_relax_to_start_coords"
	echo "Relaxing ligand structure: $cmd"
	eval $cmd

	echo "Changing chain ID of relaxed receptor and relaxed ligand back to $receptorChainId and $ligandChainId, respectively."
	cat $relaxedReceptorFile | perl -ne 'if(/^ATOM/){substr($_, 21, 1,"'$receptorChainId'")} print $_' > temp.pdb && mv temp.pdb $relaxedReceptorFile
	cat $relaxedLigandFile | perl -ne 'if(/^ATOM/){substr($_, 21, 1,"'$ligandChainId'")} print $_' > temp.pdb && mv temp.pdb $relaxedLigandFile

	echo "Concatinating both relaxed structures into a single PDB file $relaxedComplex"
	grep "^ATOM" $relaxedReceptorFile > $relaxedComplex
	echo "TER" >> $relaxedComplex
	grep "^ATOM" $relaxedLigandFile >> $relaxedComplex
    else
	echo "Concatinating receptor and ligand structures into a single PDB file $complexFile"
	grep "^ATOM" $pwd/$1 > $complexFile
	echo "TER" >> $complexFile
	grep "^ATOM" $pwd/$2 >> $complexFile

	cmd="$3 -in:file:s $complexFile -out:nstruct 1 -database $4 -docking:partners $receptorChainId""_""$ligandChainId"
	echo "Prepacking the side chain of the complex: $cmd"
	eval $cmd

	mv complex_0001.pdb $relaxedComplex
    fi	

    echo "DONE!"
    echo
    cd $pwd

fi
