#!/bin/bash

$ROSETTA3_SRC/bin/idealize.default.linuxgccrelease \
 -database $ROSETTA3_DB \
 -in::file::s 2JEL_P.pdb \
 -overwrite

mv 2JEL_P_0001.pdb 2JEL_P_idl.pdb