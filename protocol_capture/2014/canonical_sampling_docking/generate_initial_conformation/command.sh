#!/bin/bash

if [ "$MPI_RUN" == "" ]; then
   MPI_RUN=mpiexec
fi

PREFIX=""
if [ "$1" == "-n" ]; then
   PREFIX="$MPI_RUN -n $2"
fi

$PREFIX $ROSETTA_BIN/rosetta_scripts.mpi.linuxgccrelease -out:level 300 -mute all_high_mpi_rank_filebuf -database $ROSETTA_DATABASE -parser:protocol gen_initial.xml @flags_replica_dock > log

mv protAB_0001.pdb P.pdb
