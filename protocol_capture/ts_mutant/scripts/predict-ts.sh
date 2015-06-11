#!/bin/bash

usage() {
    echo "usage: $(basename $0) -protein <protein-prefix>" > /dev/stderr
    echo "analyzes Rosetta relax run output and make ts mutant predictions" > /dev/stderr
    echo "  <protein-prefix> is the name of the protein PDB file without the .pdb extension" > /dev/stderr
    exit 1
}

unset PN
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	    usage
	    ;;
	-protein)
	    PN=$2
	    shift; shift
	    ;;
	*)
	    echo "unknown argument \"$1\", aborting" > /dev/stderr
	    usage
	    ;;
    esac
done

if [ -z "$PN" ]; then
    echo "missing required argument" > /dev/stderr
    usage
fi

x=$(which java)
if [ $? -ne 0 ]; then
    echo "java is either not installed or not on your PATH" > /dev/stderr
    echo "java can be found at http://java.sun.com or installed as a Linux package" > /dev/stderr
    exit 1
fi

x=$(which psiblast)
if [ $? -ne 0 ]; then
    echo "psiblast is either not installed or not on your PATH" > /dev/stderr
    echo "NCBI BLAST+ tools can be found at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/" > /dev/stderr
    echo "Installation instructions are at http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/unix_setup.html" > /dev/stderr
    exit 1
fi

export PATH=$PATH:$(pwd)/scripts
export MODELLER_ROOT=$(pwd)
ml_src=input

echo "***"
echo "*** Creating FASTA file from PDB"
echo "***"
pdb_fasta.pl ${PN}.pdb > ${PN}.fasta

echo "***"
echo "*** Merging in extra score terms"
echo "***"
mergeall.sh ${PN}*rescore.sc

echo "***"
echo "*** Preparing score files for further processing"
echo "***"
if [ ! -e "$ml_src" ]; then mkdir $ml_src; fi
makecsv.sh -merge -wt_label non-ts ${PN}.csv > ${ml_src}/${PN}-in.csv

echo "***"
echo "*** Running PSI-BLAST (will take several minutes)"
echo "***"
if [ ! -e "pssm2" ]; then mkdir pssm2; fi
run_psiblast.sh ${PN}

echo "***"
echo "*** Preparing input file for classification"
echo "***"
prepare.sh ${ml_src}/${PN}-in.csv
mv ${ml_src}/${PN}.arff .

echo "***"
echo "*** Making predictions"
echo "***"
predict.sh ${PN}.arff
