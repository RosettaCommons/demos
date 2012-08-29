#!/usr/bin/env bash
set -e

# parameters you may wish to edit:
# num trajectories used in sequence recovery
default_n=5
n=${1:-$default_n}
# native pdbs for sequence recovery
default_nats="${my_dir}/../../inputs/min_nats/nats.list"
nats=${2:-$default_nats}

# "private" parameters you don't want to edit:
my_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
tmp_nats="inf_pro_tmp_nats.list"

for ((i=0; i<n; i++))
do
    cat "$nats"
done | sort > $tmp_nats

# break apart score.sc by protein into several score.{protein-name}.sc files
cmd="${my_dir}/score_split.py score.sc 1"
$cmd

grep total_score score.sc > top_score.sc

for i in $(ls score.*.sc)
do
    grep -v total_score $i | sort -r | head -n $n
done >> top_score.sc

# somehow this line deletes the split score files before the ls above finds them
#rm score.*.sc

${my_dir}/interfaceDesignScorer.pl -i top_score.sc > top_interface_scores

${my_dir}/topjects.sh $n

${my_dir}/seqRec.sh $tmp_nats top${n}ject.list

# Remove temporary files
rm $tmp_nats
rm top_score.sc
rm score.*.sc
