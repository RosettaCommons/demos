#!/bin/bash
# wrapper script for makecsv.awk

usage() {
    echo "usage: $(basename $0) [-merge] [-wt_label <label>] [-wt_per_group] [-d <rad>] [-ss] <input.csv>[:sasa-file.txt[:native.sc[:secondary.ss]]] ..." > /dev/stderr
    echo "usage: $(basename $0) [-merge] [-wt_label <label>] [-wt_per_group] [-ss] <input.csv>[:sasa-file.txt[:native.sc[:secondary.ss]]] ..." > /dev/stderr
    echo "wrapper script for makecsv.awk: prepares ml input (-in.csv) file from .csv name/label/mutation file(s)" > /dev/stderr
    echo "  -merge: use merge.sc/rescore.sc files" > /dev/stderr
    echo "  -wt_label: define NAT/wt line label (default is \"nc\")" > /dev/stderr
    echo "  -wt_per_group: indicates that a wt line preceeds each group (e.g. mutations at a position)" > /dev/stderr
    echo "        default behavior is to expect an initial input wt line and output it before each new group" > /dev/stderr
    echo "  -d: use local runs of radius <rad>" > /dev/stderr
    echo "  -ss: use secondary structure file" > /dev/stderr
    echo "  default sasa file is sasa-<input>.txt" > /dev/stderr
    echo "  default nat file is <input>.sc" > /dev/stderr
    echo "  assumes all .sc files are in working directory" > /dev/stderr
    exit 1
}

# defaults
wt_label="nc"
wt_per_group=0
merge=0
rad=0
unset ss

# process args
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	    usage
	    ;;
	-merge)
	    merge=1
	    shift
	    ;;
	-wt_label)
	    wt_label=$2
	    shift; shift
	    ;;
	-wt_per_group)
	    wt_per_group=1
	    shift
	    ;;
	-d)
	    rad=$2
	    shift; shift
	    ;;
	-ss)
	    ss=1
	    shift
	    ;;
	*)
	    echo "error: unknown argument $1" > /dev/stderr
	    usage
	    ;;
    esac
done

# be sure at least one input file given
if [ "$#" -eq 0 ]; then usage; fi

# process inputs
for a in $*; do
    cnt=0
    tmp=$a
    unset val
    while [ -n "$tmp" ]; do
	v="${tmp%%:*}"
	val[$cnt]=$v
	tmp="${tmp/#$v}"
	tmp="${tmp#:}"
	cnt=$(($cnt + 1))
    done
    infile=${val[0]}
    safile=${val[1]}
    natfile=${val[2]}
    ssfile=${val[3]}

    if [ -z "$infile" -o ! -f "$infile" ]; then
	echo "cannot find input file \"$infile\"" > /dev/stderr
	usage
    fi

    if [ -z "$safile" ]; then safile=sasa-${infile%.csv}.txt; fi
    if [ -z "$natfile" ]; then natfile=${infile%.csv}.sc; fi
    if [ -z "$ss" ]; then unset ssfile; else if [ -z "$ssfile" ]; then ssfile=${infile%.csv}.ss; fi; fi
    if [ -n "$rad" ]; then dstr=", d=$rad"; else dstr=""; fi

    echo "processing $a ($infile:$safile:$natfile:$ssfile)${dstr}" > /dev/stderr

    bd=$(dirname $(which $0))

    awk -f ${bd}/makecsv.awk $infile $safile $natfile ${bd}/aa.txt $wt_label $wt_per_group $merge $rad "$ssfile"
done
