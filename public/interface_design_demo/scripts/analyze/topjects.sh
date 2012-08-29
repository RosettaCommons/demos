#!/usr/bin/env bash
set -e

default=1
numJects=${1:-$default}

[ -e "toptmp" ] && rm toptmp

for i in $(ls score.*.sc) ;
do {
    tmp=${i#score.}
    prot=${tmp%.sc}
    # Check if this score file has at least $numJects number of trajectories.
    #If not, echo to stderr and exit with failure code.
    [ "$(wc -l $i | tr " " "\n" | head -n 1)" -le "$numJects" ] &&
        echo "There are not enough successful trajectories of ${prot} to" \
            "generate a list of the top $numJects for each protein." >&2 &&
        exit 1
    grep -v total_score $i | sort -r | head -n $numJects
}
done >> toptmp

awk -v pwd="$PWD" '{print pwd"/"$NF".pdb.gz"}' toptmp > "top"$numJects"ject.list"

rm toptmp
