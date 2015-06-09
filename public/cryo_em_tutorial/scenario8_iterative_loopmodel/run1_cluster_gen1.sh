#!/bin/bash

cd generation1

$ROSETTA3_SRC/bin/cluster.default.linuxgccrelease \
	-database $ROSETTA3_DB \
	-in::file::s *.pdb \
	-in:file:fullatom \
	-cluster:radius 2.5

cp c.*.0.pdb ../generation2

cd ..