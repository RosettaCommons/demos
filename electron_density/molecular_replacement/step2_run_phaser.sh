#!/bin/sh

mkdir phaser
cd phaser

phenix.phaser << eof
title 1crb
mode MR_AUTO
HKLIN ../inputs/1crb.mtz
LABIn F=FOBS_X SIGF=SIGFOBS_X
ENSEmble 1crb PDB ../templates/2qo4_mr.pdb  IDENTITY 31
COMP PROT SEQ ../inputs/1crb.fasta NUM 1
SEARCH ENSEMBLE 1crb
TOPFILES 1
ROOT 2qo4_mr.PHASER
eof
