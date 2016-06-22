#!/usr/bin/env bash

rosetta_path=/home/frank/rosetta
seq_bin=${rosetta_path}/rosetta_source/bin/sequence_recovery.linuxgccrelease
rosetta_database=${rosetta_path}/rosetta_database

my_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

nat_list=$1
redes_list=$2

$seq_bin -database $rosetta_database -native_pdb_list $nat_list -redesign_pdb_list $redes_list -parse_taskops_file ${my_dir}/seq_rec_t_ops.xml -ignore_unrecognized_res -mute all
