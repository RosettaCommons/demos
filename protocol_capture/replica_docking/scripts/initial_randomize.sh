#!/bin/bash

$ROSETTA3_BIN/rosetta_scripts.default.linuxgccrelease \
		-parser:protocol $PROTOCOL_CAPTURE/replica_docking/rosetta_inputs/randomize_infile.xml \
		-in:file:s protAB.pdb \
		-database $ROSETTA3_DB

mv protAB_0001.pdb P.pdb	


