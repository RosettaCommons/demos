#!/bin/bash
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#

if [ -z "$1" ]; then
	echo "run assembe_placed_fragment using simulated annleaing Monte Carlo sampling"
	echo "USAGE: $0 [number of jobs to run]"
	exit
fi

if [ ! -f last_dir.txt ]; then
    echo 0 > last_dir.txt
fi

last_dir=`cat last_dir.txt`
n_jobs=$1

for i in `seq $last_dir $(($last_dir+$n_jobs))`; do
	basedir=`pwd`
	if [ ! -d $i ]; then
		echo $i 
		mkdir $i
		cd $i
            ../../../rosetta/assemble_placed_fragments.static.linuxgccrelease \
            -database ../../../rosetta/rosetta_database \
            -outfile mc_sampling.out \
            -runid $i \
            -fragidx_file ../frags.idx1 \
            -densfile ../all_density.idx1 \
            -nonoverlapfile ../all_nonoverlap_scores.weighted.idx1 \
            -overlapfile ../all_overlap_scores.idx1  \
            -nmodels 500 \
            -sa_start_temp 500  \
            -sa_nsteps 200  \
            -mc_nsteps 500  \
            -null_frag_score -150  \
            -wt_dens 1.0 \
            -wt_clash 10.0 \
            -wt_closab 3.0 \
            -wt_overlap 2.0 \
            -mute all >> log &
        cd $basedir
		echo $i > last_dir.txt
    fi
done
