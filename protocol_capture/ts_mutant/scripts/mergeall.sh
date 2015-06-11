#!/bin/bash

usage() {
    echo "usage: $(basename $0) <rescore-file> [<rescore-file>...]" > /dev/stderr
    exit 1
}

if [ "$#" -eq 0 ]; then usage; fi

for x in $*; do
    echo "merging: $x"
    mergescores $x > ${x/rescore/merge}
done
