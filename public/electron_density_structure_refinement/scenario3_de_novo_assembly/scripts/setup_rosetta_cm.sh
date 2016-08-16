#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
if [ ! -f ../Step7_round2_Simulated_annealing_Monte_Carlo_sampling/average_model/average.pdb ]; then
    echo "ERROR: couldn't not find average.pdb to start with"
    exit
fi
echo Found average.pdb 
ln -s ../Step7_round2_Simulated_annealing_Monte_Carlo_sampling/average_model/average.pdb .
ln -s ../../input_files/t001_.200.9mers .
ln -s ../../input_files/t001_.200.3mers .
ln -s ../../input_files/trpv1.fasta starting.fasta
ln -s ../../input_files/transmem.mrc starting.mrc


echo Running RosettaCM now...  typically, you need to generate 500 models at least
../../rosetta/rosetta_scripts.static.linuxgccrelease @rosetta_cm_flags -mute all
