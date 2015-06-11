#!/bin/bash
# wrapper script for cutoff-sasa and gen_csv
# simplifies arguments and command line
# automatically determines aa.txt file location

function usage() {
    echo "$(basename $0) -protein <protein-prefix> -species <species> -cutoff <co>" > /dev/stderr
    echo "generates a csv file of mutations to file <protein-prefix>.pdb using all residues meeting accessibility cutoff <co>" > /dev/stderr
    echo "  -protein: specifies protein name prefix; <protein-prefix>.pdb and sasa-<protein-prefix>.txt are assumed to exist" > /dev/stderr
    echo "  -species: species name (e.g. Scer), required for later processing steps" > /dev/stderr
    echo "  -cutoff: gives maximum accessibility cutoff" > /dev/stderr
    echo "  -usage: generate this message" > /dev/stderr
    exit 1
}

unset species prefix co
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	    usage
	    ;;
	-species)
	    species=$2
	    shift; shift
	    ;;
	-protein)
	    prefix=$2
	    shift; shift
	    ;;
	-cutoff)
	    co=$2
	    shift; shift
	    ;;
	*)
	    echo "unknown argument \"$1\", aborting" > /dev/stderr
	    usage
	    ;;
    esac
done

if [ -z "$species" -o -z "$prefix" -o -z "$co" ]; then
    echo "must specify -protein, -species, and -cutoff" > /dev/stderr
    usage
fi

dn=$(dirname $0)
aafile=${dn}/aa.txt
gen=${dn}/gen_csv.awk

cutoff-sasa --list sasa-${prefix}.txt $co | awk -f $gen --species $species --pdb ${prefix}.pdb --aafile $aafile
