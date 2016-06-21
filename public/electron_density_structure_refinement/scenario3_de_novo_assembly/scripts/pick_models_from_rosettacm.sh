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
echo "cleaning rosetta_cm.out"
../../scripts/combine_silent.pl rosetta_cm.out

echo "picking models"
../../scripts/extract_lowscore_pdbs.py combined.out -p 20 -so -o 20p_combined.out --exclude_terms elec_dens_fast atom_pair_constraint -ow
../../scripts/extract_lowscore_pdbs.py 20p_combined.out -n 10 -c elec_dens_fast -l -ow

