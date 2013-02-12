#!/bin/bash

# scripts in this file is to setup target for 1bvn

cd $PROTOCOL_CAPTURE/replica_docking/start_pdbs
#RosettaDock:
	#setup target for RosettaDock as in published paper "Protein-Protein Docking with Simutaneous Optimization of Rigid-body Displacement and Side-chain Conformations", Jeffrey J. Gray et al., J. Mol. Biol. (2003)
	setup_target -method rosetta_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -disulf disulf_file -native protAB.pdb -pdb P.pdb -partners partners

#ReplicaDock:
	#setup target for ReplicaDock as described in Zhang and Lange, PLOS One 2012.
	setup_target -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -native protAB.pdb -pdb P.pdb -partners partners

	#or you can conviently copy the inputs from a previously prepared 'rosetta_dock' setup as follows:
	setup_target -method replica_dock -target udock_1bvn -target_prefix $PROTOCOL_CAPTURE/replica_docking/dock_targetlib -transfer_method rosetta_dock
