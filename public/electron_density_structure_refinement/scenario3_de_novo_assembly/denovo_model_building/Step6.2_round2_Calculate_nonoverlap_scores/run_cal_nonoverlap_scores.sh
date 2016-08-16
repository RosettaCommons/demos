#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
if [ ! -f run.lock ]; then 
    ../../../scripts/cal_NonOverlapScores.py \
        --input_dir ../../Step5_round2_Place_fragments_into_density/combined_round1_round2_frags/ 


    if [ -f run.lock ]; then # this means the previous step runs
        ../../../scripts/clean_scorefile.py \
            --selected_frags_path  ../../Step5_round2_Place_fragments_into_density/combined_round1_round2_frags/ \
            -w \
            -t nonoverlap \
            -nw 10  
    fi
fi
