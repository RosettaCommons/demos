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
# set up input files for SAMC runs
rt=$PWD

cd ../Step5_round2_Place_fragments_into_density/
    echo get frags.idx1, all_density.idx1
    ../../scripts/get_idx1.py -p combined_round1_round2_frags/ -c combined_round1_round2_frags.sc
cd $rt

cd ../Step6.1_round2_Calculate_overlap_scores 
    echo get all_overlap_scores.idx1
    cat ./*/*idx > all_overlap_scores.idx1
cd $rt

cd ../Step6.2_round2_Calculate_nonoverlap_scores
    echo get all_nonoverlap_scores.weighted.idx1
    cat ./*/*idx > all_nonoverlap_scores.weighted.idx1
cd $rt

# link all the idx1 files here
ln -s ../*round2*/*idx1 .

