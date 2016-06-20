#!/bin/bash
#
# creates intermediate .csv and final .arff file for a single command line arg

usage() {
    echo "usage: $(basename $0) [-v <n>] [-c <n>] [-prefix <prefix>] <file-in.csv> [extra-args]" > /dev/stderr
    echo "creates <file>.csv and <file>.arff" > /dev/stderr
    echo "  -v: version number (1-3)" > /dev/stderr
    echo "  -c: number of output classes (2 or 3)" > /dev/stderr
    echo "  -prefix: give an explicit prefix (default is to strip \"-in.csv\" from input file name" > /dev/stderr
    echo "  extra-args: passed directly to jp.sh" > /dev/stderr
    exit 1
}

prep_ver=2
prep_cat=2
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-v)
	    prep_ver=$2
	    if [ "$prep_ver" -lt 1 -o "$prep_ver" -gt 3 ]; then
		echo "error: invalid -v argument \"$prep_ver\"" > /dev/stderr
		usage
	    fi
	    shift; shift
	    ;;
	-c)
	    prep_cat=$2
	    if [ "$prep_cat" -lt 2 -o "$prep_cat" -gt 3 ]; then
		echo "error: invalid -c argument \"$prep_cat\"" > /dev/stderr
		usage
	    fi
	    shift; shift
	    ;;
	-prefix)
	    prefix=$2
	    shift; shift
	    ;;
	-usage)
	    usage
	    ;;
	*)
	    echo "error: unrecognized argument $1" > /dev/stderr
	    usage
	    ;;
    esac
done

infile=$1
extras=${*:2}  # args 2:end

if [ ! -f "$infile" ]; then
    echo "error: could not find input file $infile" > /dev/stderr
    usage
fi

if [ -z "$prefix" ]; then prefix="${infile%-in.csv}"; fi

if [ $prep_cat -eq 2 ]; then
    labels="non-ts,ts"
else
    labels="nc,ts,lof"
fi

outfile=${prefix}.csv
#arff_opt="-v${prep_ver}c${prep_cat}"
if [ $prep_ver -eq 1 ]; then jp_args="-noconvertlabels"; else unset jp_args; fi
jp.sh -v${prep_ver} $jp_args $extras $infile > $outfile
convert-arff-new.sh $outfile $labels
# add source csv for later indexing of predictions
sed -i -e '1i\
% from '"$infile" ${outfile/.csv/.arff}
