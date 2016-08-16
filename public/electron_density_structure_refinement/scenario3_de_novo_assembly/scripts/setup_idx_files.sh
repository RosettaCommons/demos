#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
# set up input files for SAMC runs
rt=$PWD

cd ../Step1_Place_fragments_into_density 
    echo get candidate_fragment_placements.sc 
    cd condor_jobs/
        cat ./*/clustering_selected.sc | grep -v "^#" > candidate_fragment_placements.sc 
    cd ../
    echo get frags.idx1, all_density.idx1
    ../../scripts/get_idx1.py -p candidate_fragment_placements/ -c condor_jobs/candidate_fragment_placements.sc
cd $rt

cd ../Step2.1_Calculate_overlap_scores 
    echo get all_overlap_scores.idx1
    cat ./*/*idx > all_overlap_scores.idx1
cd $rt

cd ../Step2.2_Calculate_nonoverlap_scores
    echo get all_nonoverlap_scores.weighted.idx1
    cat ./*/*idx > all_nonoverlap_scores.weighted.idx1
cd $rt

# link all the idx1 files here
ln -s ../*/*idx1 .

