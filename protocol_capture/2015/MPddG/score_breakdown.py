#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file: compute_ddG.py
##
## @brief:   Compute ddGs of mutation
## @details: Use the Rosetta membrane framework to compute the ddG of unfolding of 
## a membrane protein in Rosetta (uses packer, mutate.py from Evan Baugh)
##
## @author: Rebecca F. Alford (rfalford12@gmail.com)
## @author: JKLeman (julia.koehler1982@gmail.com)

# Tools
import sys, os
import commands
import random
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

# Rosetta-specific imports
import rosetta.protocols.membrane
from rosetta import Pose
from rosetta import create_score_function
from rosetta import TaskFactory
from rosetta.utility import vector1_bool
from rosetta import aa_from_oneletter_code
from rosetta import PackRotamersMover
from rosetta.core.pose import PDBInfo
from rosetta.core.chemical import VariantType

###############################################################################

## @brief Main - Add Membrane to Pose, Compute ddG
def main( args ):
	
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    #input options
    parser.add_option('--in_pdb', '-p',
       action="store",
       help="Input PDB file.", )

    parser.add_option('--in_span', '-s',
       action="store",
       help="Input spanfile.", )
				  
    parser.add_option('--out', '-o',
       action="store", default='ddG.out',
       help="Output filename with pose residue numbering. Default: 'ddG.out'", )
									
    parser.add_option('--res', '-r',
       action="store",
       help="Pose residue number to mutate.", )

    parser.add_option('--mut', '-m',
       action="store",
       help="One-letter code of residue identity of the mutant. Example: A181F would be 'F'", )

    parser.add_option('--repack_radius', '-a', 
        action="store", default=0, 
        help="Repack the residues within this radius",)

    parser.add_option('--output_breakdown', '-b', 
        action="store", default="scores.sc", 
        help="Output mutant and native score breakdown by weighted energy term into a scorefile", )

    parser.add_option('--include_pH', '-t', 
        action="store", default=0,
        help="Include pH energy terms: pH_energy and fa_elec. Default false.", )

    parser.add_option('--pH_value', '-q', 
        action="store", default=7,
        help="Predict ddG and specified pH value. Default 7. Will not work if include pH is not passed", )

    #parse options
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check the required inputs (PDB file, spanfile) are present
    if ( not Options.in_pdb or not Options.in_span or not Options.res ):
	    sys.exit( "Must provide flags '-in_pdb', '-in_span', and '-res'! Exiting..." )

    # Initialize Rosetta options from user options. Enable pH mode if applicable
    rosetta_options = ""
    standard_options = "-mp:setup:spanfiles " + Options.in_span +  " -run:constant_seed -in:ignore_unrecognized_res"
    if ( Options.include_pH ): 
        print Options.pH_value
        if ( float( Options.pH_value ) < 0 or float(Options.pH_value) > 14 ): 
            sys.exit( "Specified pH value must be between 0-14: Exiting..." )
        else: 
            pH_options = " -pH_mode -value_pH " + str(Options.pH_value)
            rosetta_options = standard_options + pH_options
    else: 
        rosetta_options = standard_options

    # Initialize Rosetta based on user inputs
    rosetta.init( extra_options=rosetta_options )
	
    # Load Pose, & turn on the membrane
    pose = pose_from_pdb( Options.in_pdb )

    # Add Membrane to Pose
    add_memb = rosetta.protocols.membrane.AddMembraneMover()
    add_memb.apply( pose )
    
    # Setup in a topology based membrane
    init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
    init_mem_pos.apply( pose )

    # check the user has specified a reasonable value for the pH


    # Create a membrane energy function enabled by pH mode
    # Includes two terms not standard in the smoothed energy function: pH energy
    # and fa_elec
    sfxn = rosetta.core.scoring.ScoreFunction()
    sfxn = create_score_function( "mpframework_pHmode_fa_2015")

    print compute_ddG( native, mutant, sfxn )

###############################################################################

## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, mutated_pose, sfxn ): 

    # Score Native Pose
    native_score = sfxn( pose )

    # Score Mutated Pose
    mutant_score = sfxn( mutated_pose )

	# return scores
    return aa, round( mutant_score, 3 ), round( native_score, 3 ), round ( mutant_score - native_score, 3 )

if __name__ == "__main__" : main(sys.argv)

