#!/bin/bash
if [ $# -ne 8 ]
then
    echo
    echo "Demo execution of cross-link guided protein docking on the example of the IgBP1-PP2AA complex. In total the paths to 8 applications and libraries must be provided."
    echo
    echo "Usage: $0 1-path2rosettaBinDirectory 2-path2rosettaDatabase 3-path2XwalkBin 4-path2naccess 5-path2cX.jar 6-path2colt.jar 7-path2cdk-1.0.2.jar 8-path2R"
    echo
    echo "Example: $0 rosetta/3.4/rosetta_source/bin/ rosetta/3.4/rosetta_database/ /src/java/xwalk_v0.5/bin/ /bin/naccess2.1.1/naccess /src/java/cX.jar /src/java/colt.jar /src/java/cdk-1.0.2.jar /usr/bin/R"
    echo
else
    rBinDir=$1
    rosettaDB=$2
    xwalkApp=$3
    naccessApp=$4
    cXlib=$5
    coltLib=$6
    cdkLib=$7
    rApp=$8

    # check that all applications and libraries really exist
    exit=0
    if [ ! -d $rBinDir ]; then echo "    ERROR: ROSETTA bin directory $rBinDir does not exist. Please verify that it exists."; exit=1; fi 
    if [ ! -d $rosettaDB ]; then echo "    ERROR: ROSETTA database directory $rosettaDB does not exist. Please verify that it exists."; exit=1; fi 
    if [ ! -d $xwalkApp ]; then echo "    ERROR: Xwalk bin directory does not exist at $xwalkApp. Please verify that it exists."; exit=1; fi
    if [ ! -f "$xwalkApp/Xwalk.class" ]; then echo "    ERROR: Xwalk exectuable does not exist in the directory $xwalkApp. Please verify that it exists."; exit=1; fi
    if [ ! -f $naccessApp ]; then echo "    ERROR: NACCESS executable does not exist at $naccessApp. Please verify that it exists."; exit=1; fi
    if [ ! -f $cXlib ]; then echo "    ERROR: CleftXplorer JAVA library file cX.jar does not exist at $cXlib. Please verify that it exists."; exit=1; fi
    if [ ! -f $coltLib ]; then echo "    ERROR: colt.jar JAVA library file does not exist at $coltLib. Please verify that it exists."; exit=1; fi
    if [ ! -f $cdkLib ]; then echo "    ERROR: cdk-1.0.2.jar JAVA library file does not exist at $cdkLib. Please verify that it exists."; exit=1; fi
    if [ ! -f $rApp ]; then if [[ `type R` ]] ; then rApp="R"; else echo "    ERROR: R executable does not exist at $rApp. Please verify that it exists."; exit=1; fi; fi
    if [[ ! `type java` ]]; then echo "    ERROR: JAVA could not be found on your operating system. Please download and install it from http://java.com"; exit=1; fi
    if [ $exit -eq 1 ]; then exit; fi

    # get the absolute paths for all applications and libraries.
    pwd=`pwd`
    cd $rBinDir; rBinDir=`pwd`; cd $pwd
    cd $rosettaDB; rosettaDB=`pwd`; cd $pwd
    cd $xwalkApp; xwalkApp=`pwd`; cd $pwd
    cd `dirname $naccessApp`; naccessApp=`pwd`"/"`basename $naccessApp`; cd $pwd
    cd `dirname $cXlib`; cXlib=`pwd`"/"`basename $cXlib`; cd $pwd
    cd `dirname $coltLib`; coltLib=`pwd`"/"`basename $coltLib`; cd $pwd
    cd `dirname $cdkLib`; cdkLib=`pwd`"/"`basename $cdkLib`; cd $pwd
    # set absolute path to R only if user has given a correct full path of R, otherwise use system wide R installation.
    if [[ ! `type R` ]]; then cd `dirname $rApp`; rApp=`pwd`"/"`basename $rApp`; cd $pwd; fi

    # find prepack and docking protocol within the ROSETTA bin directory. 
    prepackApp=`ls $rBinDir/docking_prepack_protocol.* | head -1`
    dockingApp=`ls $rBinDir/docking_protocol.* | head -1`
    if [[ ! -f $prepackApp ]]; then echo "    ERROR: docking_prepack_protocol does not exist in ROSETTA bin directory $rBinDir. Please verify that it exists."; exit=1; fi
    if [[ ! -f $dockingApp ]]; then echo "    ERROR: docking_protocol does not exist in ROSETTA bin directory $rBinDir. Please verify that it exists."; exit=1; fi
    if [ $exit -eq 1 ]; then exit; fi


    # the relative path to the receptor PDB file. Should be the larger protein among the binding partners.
    receptorFile="inputs/igbp1.pdb"

    # the relative path to the ligand PDB file. Should be the smaller protein among the binding partners.
    ligandFile="inputs/pp2aa.pdb"

    # native protein complex, required for RMSD calculation. If no native complex exists, use a concatinated PDB file of the receptor and ligand coordinates.
    nativeFile="inputs/igbp1-pp2aa-bestModel-complex.pdb"

    # list of inter-protein cross-links in ROSETTA constraint file format.
    cstFile="inputs/igbp1-pp2aa-xls.cst"

    # list of intra-protein cross-links, inter-protein cross-links and mono-links in Xwalk distance file format.
    xlsFile="inputs/igbp1-pp2aa-xls.txt"

    # do relax instead of side chain prepacking
    doRelax=0

    # number of global docking models
    globalN=5

    # number of local docking models per cluster representative
    localN=5

    # remove models with too small interfaces
    useInterfaceFilter=0

    echo "--------------"
    echo "10. PREPACKING"
    echo "--------------"
    ./scripts/1-prepack_or_relax.sh $receptorFile $ligandFile $prepackApp $rosettaDB $doRelax
    echo "-----------------"
    echo "9. GLOBAL DOCKING"
    echo "-----------------"
    ./scripts/2-globaldock.sh $dockingApp $rosettaDB $globalN $nativeFile $cstFile 
    echo "-------------------------------"
    echo "8. XWALK DISTANCE CALCULATION I"
    echo "-------------------------------"
    ./scripts/3-xwalk.sh $xwalkApp output/docking/global output/cluster/global $xlsFile 
    echo "-----------------------------"
    echo "7. ASSESSING INTERFACE SIZE I"
    echo "-----------------------------"
    ./scripts/4-naccess.sh $naccessApp $useInterfaceFilter output/cluster/global 
    echo "----------------"
    echo "6. QT-CLUSTERING"
    echo "----------------"
    ./scripts/5-qtcluster.sh $cXlib $coltLib $cdkLib 
    echo "---------------------------"
    echo "5. LOCAL REFINEMENT DOCKING"
    echo "---------------------------"
    ./scripts/6-localdock.sh $dockingApp $rosettaDB $localN $nativeFile $cstFile 
    echo "--------------------------------"
    echo "4. XWALK DISTANCE CALCULATION II"
    echo "--------------------------------"
    ./scripts/7-xwalk.sh $xwalkApp output/docking/local output/cluster/local $xlsFile 
    echo "------------------------------"
    echo "3. ASSESSING INTERFACE SIZE II"
    echo "------------------------------"
    ./scripts/8-naccess.sh $naccessApp $useInterfaceFilter output/cluster/local 
    echo "------------------"
    echo "2. RMSD-CLUSTERING"
    echo "------------------"
    ./scripts/9-rmsdcluster.sh $cXlib $coltLib $cdkLib $rApp
    echo "-----------------------"
    echo "1. PREDICTING INTERFACE"
    echo "-----------------------"
    ./scripts/10-interface.sh $xwalkApp $naccessApp
    echo "-------"
    echo "0. DONE"
    echo "-------"
fi
