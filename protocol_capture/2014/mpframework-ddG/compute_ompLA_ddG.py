#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file: compute_ompLA_ddG.py
##
## @brief: Compute ddG of mutation in ompLA 
## @details: Use the Rosetta membrane framework to compute the ddG of unfolding of 
## a membrane protein in Rosetta (uses packer, mutate.py from Evan Baugh) and compare
## to experimentally meausred values in Moon & Fleming, 2011
##
## Citation: Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from 
## transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. Available: 
## http://www.pnas.org/content/early/2011/05/16/1103979108.
##
## @author: Rebecca F. Alford (rfalford12@gmail.com)

from rosetta import *

# tools
import ddG
import sys, os
import commands

# rosetta
import rosetta.protocols.membrane

## Default Values
repack_radius = 0.0

## @brief Main - Add Membrane to Pose, Compute ddG
def main( argv ):

    rosetta.init( extra_options="-membrane_new:setup:spanfiles inputs/1qd6_tr.span -run:constant_seed -in:ignore_unrecognized_res" )

    # Load Pose, & turn on the membrane
    pose = pose_from_pdb( "inputs/1qd6.pdb" );
    sfxn = create_score_function( "fa_menv_smooth_2014" );
    
    # Add Membrane to Pose
    add_memb = rosetta.protocols.membrane.AddMembraneMover()
    add_memb.apply( pose )

    # Setup in a topology based membrane
    init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
    init_mem_pos.apply( pose )

    # Compute mutations 
    AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    ddGs = []
    for aa in AAs:
        print compute_ddG( pose, sfxn, 181, aa  )

if __name__ == "__main__" : main(sys.argv)

