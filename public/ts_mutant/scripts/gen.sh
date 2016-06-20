#!/bin/bash

usage() {
    echo "usage: $(basename $0) -native native.pdb[.gz] [-bin <mini_bin_path>] [-db <mini_db_path>] input.pdb[.gz] ..." > /dev/stderr
    exit 1
}

unset nat
mini_bin=$MINI_BIN
mini_db=$MINI_DB
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	    usage
	    ;;
	-native)
	    nat=$2
	    shift; shift
	    ;;
	-bin)
	    mini_bin=$2
	    shift; shift
	    ;;
	-db)
	    mini_db=$2
	    shift; shift
	    ;;
	*)
	    echo "unknown argument \"$1\", aborting" > /dev/stderr
	    usage
	    ;;
	esac
done

if [ -z "$nat" ]; then
    echo "must specify -native" > /dev/stderr
    usage
fi
if [ -z "$mini_bin" ]; then
    echo "must specify -bin or define \$MINI_BIN" > /dev/stderr
    usage
fi
if [ -z "$mini_db" ]; then
    echo "must specify -db or define \$MINI_DB" > /dev/stderr
    usage
fi

tmp=${nat%.gz}
prefix_in=${tmp%.pdb}

for fn in $nat $*; do
    bnt=${fn%.gz}        # remove .gz if present
    bn=${bnt%.pdb}       # remove .pdb in all cases
    if [ -z "$prefix_in" ]; then
	prefix=${bn%-*}      # remove -.RESI.
    else
	prefix=$prefix_in
    fi
    id="${bn#$prefix-}"  # isolate .RESI. (e.g. V33A or 101)
    if [ $bn == $prefix ]; then  # there was no .RESI.
	m=WT
    else
	m=$id
    fi
    of=$prefix-${m}.sh
    if [ "${fn%.gz}" != "$fn" ]; then nmex=".pdb_????.pdb.gz"; else nmex="_????.pdb.gz"; fi
#    echo "fn: $fn, bnt: $bnt, bn: $bn, prefix: $prefix, id: $id, m: $m, of: $of"
    echo "creating $of ($m) from $fn"
    echo "$mini_bin/relax.linuxgccrelease -database $mini_db -s $fn -native $nat -nstruct 50 -relax:fast -out:file:scorefile $prefix-${m}.sc -out:pdb_gz" > $of
    echo "$mini_bin/score.linuxgccrelease -database $mini_db -s ${bn}${nmex} -in:file:native $nat -in:file:fullatom -out:file:scorefile $prefix-${m}rescore.sc" >> $of
    chmod 755 $of
    # generate native score script
    if [ $m == WT ]; then
	sf=${prefix}.sh
	echo "creating native score script $sf"
	echo "$mini_bin/score.linuxgccrelease -database $mini_db -s $nat -out:file:scorefile ${prefix}.sc" > $sf
	chmod 755 $sf
    fi
    shift
done
