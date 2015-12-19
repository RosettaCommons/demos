#!/bin/bash

# scripts in this file are to generate the run directories with CS-Rosetta Automated setup toolbox

# ReplicaDock, using temperatures [2.0 3.0 5.0]
setup_run -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock -job interactive -extras mpi -score interchain_cen -nstruct 1 -protocol rep_cen -xml uniform -n_replica 3

# ReplicaDock-LoT, using temperatures [0.6 0.8 1.0 1.2 1.5 2.0 2.5] + min_score
setup_run -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/replica_dock_LoT -job interactive -extras mpi -score interchain_cen -min_score -36.519 -nstruct 1 -protocol rep_cen -xml uniform -n_replica 7

# RosettaDock, only low resolution phase
setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/rosetta_dock -job interactive -extras mpi -protocol centroid -batches 2 -score interchain_cen -nstruct 25

# refinement of the low-resolution ensembles produced by ReplicaDock
setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_replica_dock -job interactive -extras mpi -protocol refine -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 -start $PROTOCOL_CAPTURE/replica_docking/example_runs/replica_dock/udock_1bvn/run

# refinement of the low-resolution ensembles produced by RosettaDock
setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/refine_rosetta_dock -job interactive -extras mpi -protocol refine -pattern "low_decoys_*out" -prefix refine -score docking -nstruct 1 -start $PROTOCOL_CAPTURE/replica_docking/example_runs/rosetta_dock/udock_1bvn/run

# generate RelaxedNative ensembles with 1000 decoys from the superimposed native structure protAB.pdb, used as reference in final analysis
setup_run -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -dir $PROTOCOL_CAPTURE/replica_docking/test_runs/relax_native -job interactive -extras mpi -protocol refine -out relax_native.out -score docking -nstruct 1000
