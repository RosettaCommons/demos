#!/bin/bash

$ROSETTA3_SRC/bin/minirosetta.default.linuxgccrelease \
	-run:protocol broker \
	-database $ROSETTA3_DB/ \
	-broker:setup setup_dens.tpb \
	-database /work/dimaio/minirosetta_database \
	-skip_stages 1 2 \
	-seq_sep_stages 1 1 1 \
	-short_frag_cycles 1 \
	-scored_frag_cycles 1 \
	-increase_cycles 0.1 \
	-ramp_chainbreaks \
	-sep_switch_accelerate 0.8 \
	-skip_convergence_check \
	-overlap_chainbreak \
	-fail_on_bad_hbond false \
	-edensity::mapfile 1nsf.5A.mrc \
	-edensity::mapreso 8.0 \
	-edensity::grid_spacing 4.0 \
	-abinitio::stage3a_patch ./phase3.patch \
	-abinitio::stage3b_patch ./phase3.patch \
	-abinitio::stage4_patch  ./phase4.patch \
	-out:file:scorefile score.fsc \
	-nstruct 1 -overwrite
