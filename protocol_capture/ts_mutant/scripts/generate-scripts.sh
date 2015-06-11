#!/bin/bash

usage() {
    echo "usage: $(basename $0) -protein <protein-prefix> -species <species-abbrev> -cutoff <cutoff>" > /dev/stderr
    echo "       -mini_bin <mini_bin_path> -mini_db <mini_db_path>" > /dev/stderr
    echo "prepares shell scripts to execute Rosetta runs for ts mutant prediction" > /dev/stderr
    echo "  <protein-prefix> is the name of the protein PDB file without the .pdb extension" > /dev/stderr
    echo "  <species-abbrev> is the species short name, e.g. Scer" > /dev/stderr
    echo "  <cutoff> is the ACC cutoff for native residues to mutate, typically 10" > /dev/stderr
    echo "  <mini_bin_path> is the path of the directory containing Rosetta binaries *on the machine" > /dev/stderr
    echo "      where Rosetta will be run*" > /dev/stderr
    echo "  <mini_db_path> is the path of the Rosetta database *on the machine where Rosetta will be run*" > /dev/stderr
    exit 1
}

unset PN species cutoff mini_bin mini_db
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	    usage
	    ;;
	-protein)
	    PN=$2
	    shift; shift
	    ;;
	-species)
	    species=$2
	    shift; shift
	    ;;
	-cutoff)
	    cutoff=$2
	    shift; shift
	    ;;
	-mini_bin)
	    mini_bin=$2
	    shift; shift
	    ;;
	-mini_db)
	    mini_db=$2
	    shift; shift
	    ;;
	*)
	    echo "unknown argument \"$1\", aborting" > /dev/stderr
	    usage
	    ;;
    esac
done

if [ -z "$PN" -o -z "$species" -o -z "$cutoff" -o -z "$mini_bin" -o -z "$mini_db" ]; then
    echo "missing required argument" > /dev/stderr
    usage
fi

x=$(which probe)
if [ $? -ne 0 ]; then
    echo "probe is either not installed or not on your PATH" > /dev/stderr
    echo "probe can be found at http://kinemage.biochem.duke.edu/software/probe.php" > /dev/stderr
    exit 1
fi

x=$(which pymol)
if [ $? -ne 0 ]; then
    echo "pymol is either not installed or not on your PATH" > /dev/stderr
    echo "pymol can be found at http://www.pymol.org" > /dev/stderr
    exit 1
fi

export PATH=$PATH:$(pwd)/scripts

echo "***"
echo "*** Generating surface area accessibility file"
echo "***"
gen-acc ${PN}.pdb > sasa-${PN}.txt

echo "***"
echo "*** Generating master list of mutations"
echo "***"
gen_csv.sh -protein $PN -species $species -cutoff $cutoff

echo "***"
echo "*** Generating pdb file for each mutation"
echo "***"
mkmut-csv ${PN}.csv

echo "***"
echo "*** Generating relax run script for each mutation plus native"
echo "***"
gen.sh -bin $mini_bin -db $mini_db -native ${PN}.pdb ${PN}-*.pdb
