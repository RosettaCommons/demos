#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) All the files in this directory and sub-directories are part of the Rosetta software
# (c) suite and are made available under license.  The Rosetta software is developed by the
# (c) contributing members of the Rosetta Commons. For more information, see
# (c) http://www.rosettacommons.org. Questions about this can be addressed to University of
# (c) Washington UW TechTransfer, email: license@u.washington.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
if [ ! -f run.lock ]; then 
    ../../../scripts/cal_OverlapScores.py \
        --input_dir ../../Step5_round2_Place_fragments_into_density/combined_round1_round2_frags/ 


    if [ -f run.lock ]; then # this means the previous step runs
        ../../../scripts/clean_scorefile.py \
            --selected_frags_path  ../../Step5_round2_Place_fragments_into_density/combined_round1_round2_frags/ \
            -t overlap 
    fi
fi
