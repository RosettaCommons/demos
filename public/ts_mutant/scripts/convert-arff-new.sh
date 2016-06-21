#!/bin/bash
#
# 

usage() {
    echo "usage: $(basename $0) <infile.csv> <labels>" > /dev/stderr
    echo "makes an ARFF file from the given CSV file" > /dev/stderr
    echo "parses first line of CSV file to create ARFF header" > /dev/stderr
    echo "  <infile.csv>: file to parse for header (product of -in.csv file)" > /dev/stderr
    echo "  <labels>: comma-separated labels for final label field" > /dev/stderr
    exit 1
}

if [[ "$1" =~ "--*usage" ]]; then usage; fi
if [ "$#" -ne 2 ]; then echo "error: arguments incorrect" > /dev/stderr; usage; fi
if [ ! -f "$1" ]; then echo "error: cannot find file \"$1\"" > /dev/stderr; usage; fi

infile=$1
labels=$2
outfile=${infile/.csv/.arff}

echo "creating $outfile"

(
    echo "@relation $(basename ${infile%.csv})"
    echo ""

    head -1 $infile | sed 's/,/\n/g' | while read field; do
	case $field in
	    RES|MUT|PT|SPECIES)
	    echo "@attribute $field string"
	    ;;
	    PT_val)
	    echo "@attribute $field {non-ts,ts}"
	    ;;
	    *)
	    echo "@attribute $field numeric"
	    ;;
	esac
    done

    echo ""
    echo "@data"
    sed -e '1d' $infile
) > $outfile
