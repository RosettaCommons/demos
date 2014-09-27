#!/usr/bin/python
# :noTabs=true:  
# (c) Copyright Rosetta Commons Member Institutions. 
# (c) This file is part of the Rosetta software suite and is made available 
# (c) under license. 
# (c) The Rosetta software is developed by the contributing members of the 
# (c) Rosetta Commons. 
# (c) For more information, see http://www.rosettacommons.org. 
# (c) Questions about this can be addressed to University of Washington UW 
# (c) TechTransfer, email: license@u.washington.edu.

## @file: membrane_ddG.py
##
## @brief: Compute ddG of mutation for a membrane protein
## @details: Use the Rosetta membrane framework to compute the ddG of unfolding of 
## a membrane protein in Rosetta (uses packer, mutate.py from Evan Baugh)
##
## @author: Rebecca F. Alford (rfalford12@gmail.com)

from rosetta import *
rosetta.init()

# tools
import mutate
import sys, os
import commands

# rosetta
import protocols.membrane

## Default Values
repack_radius = 0.0

## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, sfxn, resnum, aa ): 

	# Score Native Pose
	native_score = sfxn( pose )

	# Perform Mutation at residue <resnum> to amino accid <aa>
	mutated_pose = mutate_residue( pose, resnum, aa, repack_radius, sfxn )

	# Score Mutated Pose
	mutant_score = sfxn( mutated_pose )

	# Compute ddG of mutation
	ddG = mutant_score - native_score

	# Print result
	print ddG

## @brief Main - Add Membrane to Pose, Compute ddG
def main( argv ):

	init( extra_options="-membrane_new:setup:spanfile -run:constant_seed" )

	# load pose
	pose = pose_from_pdb("1py6_tr_native.pdb");
	AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

	for aa in AAs:
		calculateDDG( pose, 181, aa, "score12")	

if __name__ == "__main__" : main(sys.argv)

