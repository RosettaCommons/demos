#!/bin/bash

	$ROSETTA3_BIN/score_jd2.default.linuxgccrelease -in:file:s protAB.pdb > log
	cat log | grep "disulfide between" | awk '{print $6,$7}' > disulf_file	
