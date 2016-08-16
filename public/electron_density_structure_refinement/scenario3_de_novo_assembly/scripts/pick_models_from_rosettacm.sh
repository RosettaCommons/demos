#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
echo "cleaning rosetta_cm.out"
../../scripts/combine_silent.pl rosetta_cm.out

echo "picking models"
../../scripts/extract_lowscore_pdbs.py combined.out -p 20 -so -o 20p_combined.out --exclude_terms elec_dens_fast atom_pair_constraint -ow
../../scripts/extract_lowscore_pdbs.py 20p_combined.out -n 10 -c elec_dens_fast -l -ow

