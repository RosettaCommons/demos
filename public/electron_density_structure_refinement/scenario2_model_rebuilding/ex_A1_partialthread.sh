#!/bin/bash

~/rosetta_workshop/rosetta/main/source/bin/partial_thread.linuxgccrelease \
 -database ~/rosetta_workshop/rosetta/main/database/ \
 -in::file::fasta t20s.fasta \
 -in::file::alignment 20S_1iru.ali \
 -in::file::template_pdb 1iruAH_aln.pdb
