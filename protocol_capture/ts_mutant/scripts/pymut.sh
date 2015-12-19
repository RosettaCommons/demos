#!/bin/bash
# front end to mk.py
#   generates mutations of residue at a position <resi> in <file.pdb>
#   passes all args to mk.py; current behavior is to generate named mutations, or all if none given

usage() {
    echo "usage: $(basename $0) <file.pdb> <resi> [<resn> ...]" > /dev/stderr;
    echo "       <resi> can optionally start with a one-char native res description" > /dev/stderr;
    echo "       <resn> is one or more one-char mutation description" > /dev/stderr;
    echo "       default is all residues if no <resn> arguments given" > /dev/stderr;
    exit 1
}

if [[ "$#" -lt 2 || "$1" =~ "-+usage" ]]; then usage; fi

dn=$(dirname $0)
pymol -qcr ${dn}/mk.py -- $*
