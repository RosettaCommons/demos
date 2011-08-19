#!/bin/bash

$ROSETTA3_SRC/bin/idealize.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s S_0001.pdb \
 -overwrite

mv S_0001.pdb S_0001_idl.pdb