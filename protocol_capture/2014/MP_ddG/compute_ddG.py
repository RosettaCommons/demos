#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file: compute_ddG.py
##
## @brief: 	 Compute ddGs of mutation
## @details: Use the Rosetta membrane framework to compute the ddG of unfolding of 
## a membrane protein in Rosetta (uses packer, mutate.py from Evan Baugh)
##
## @author: Rebecca F. Alford (rfalford12@gmail.com)
## @author: JKLeman (julia.koehler1982@gmail.com)


# tools
import sys, os
import commands
import random
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

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

    #parse options
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # check whether all inputs are there
    if ( not Options.in_pdb or not Options.in_span or not Options.res ):
	    sys.exit( "Must provide flags '-in_pdb', '-in_span', and '-res'! Exiting..." )

    # importing Rosetta
    import rosetta.protocols.membrane
    from rosetta import Pose
    from rosetta import create_score_function
    from rosetta import TaskFactory
    from rosetta.utility import vector1_bool
    global vector1_bool
    from rosetta import aa_from_oneletter_code
    from rosetta import PackRotamersMover
    from rosetta.core.pose import PDBInfo

    ## Default Values
    global repack_radius
    repack_radius = 0.0

    # initialize Rosetta
    rosetta.init( extra_options="-membrane_new:setup:spanfiles " + Options.in_span + " -run:constant_seed -in:ignore_unrecognized_res" )
	
    # Load Pose, & turn on the membrane
    pose = pose_from_pdb( Options.in_pdb );
    sfxn = create_score_function( "mpframework_smooth_fa_2014" );
#    sfxn = create_score_function( "mpframework_pHmode_fa_2014" );
    
    # Add Membrane to Pose
    add_memb = rosetta.protocols.membrane.AddMembraneMover()
    add_memb.apply( pose )
	
    # Setup in a topology based membrane
    init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
    init_mem_pos.apply( pose )
	
    # Compute mutations
    if ( Options.mut ):
        with file( Options.out, 'a' ) as f:
            ddGs = compute_ddG( pose, sfxn, int( Options.res ), Options.mut )
            f.write( Options.in_pdb + " " + Options.res + " " + str(ddGs[0]) + " " + str(ddGs[1]) + " " + str(ddGs[2]) + " " + str(ddGs[3]) + "\n" )
	    f.close
    else:
        #AAs = ['D', 'E']
        AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        for aa in AAs:
            with file( Options.out, 'a' ) as f:
                ddGs = compute_ddG( pose, sfxn, int( Options.res ), aa )
                f.write( str(ddGs[0]) + " " + str(ddGs[1]) + " " + str(ddGs[2]) + " " + str(ddGs[3]) + "\n" )
            f.close

###############################################################################

## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, sfxn, resnum, aa ): 

    # Score Native Pose
    native_score = sfxn( pose )

    # Perform Mutation at residue <resnum> to amino acid <aa>
    mutated_pose = mutate_residue( pose, resnum, aa, repack_radius, sfxn )

    # Adjust for pH
    #adjust_for_pH( pose, resnum )

    # Score Mutated Pose
    mutant_score = sfxn( mutated_pose )

	# return scores
    return aa, round( mutant_score, 3 ), round( native_score, 3 ), round ( mutant_score - native_score, 3 )

###############################################################################

# @brief Replace the residue at <resid> in <pose> with <new_res> and allows
# repacking within a given <pack_radius> 
def mutate_residue( pose , mutant_position , mutant_aa ,
        pack_radius = 0.0 , pack_scorefxn = '' ):

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign( pose )

    # Create a talaris sfxn by default
    if not pack_scorefxn:
        pack_scorefxn = create_score_function( 'fa_menv_smooth_2014' )

    # Create a packer task (standard)
    task = TaskFactory.create_packer_task( test_pose )


    # the Vector1 of booleans (a specific object) is needed for specifying the
    #    mutation, this demonstrates another more direct method of setting
    #    PackerTask options for design
    aa_bool = vector1_bool()

    # PyRosetta uses several ways of tracking amino acids (ResidueTypes)
    # the numbers 1-20 correspond individually to the 20 proteogenic amino acids
    # aa_from_oneletter returns the integer representation of an amino acid
    #    from its one letter code
    # convert mutant_aa to its integer representation
    mutant_aa = aa_from_oneletter_code( mutant_aa )

    # mutation is performed by using a PackerTask with only the mutant
    #    amino acid available during design
    # to do this, construct a Vector1 of booleans indicating which amino acid
    #    (by its numerical designation, see above) to allow
    for i in range( 1 , 21 ):
        # in Python, logical expression are evaluated with priority, thus the
        #    line below appends to aa_bool the truth (True or False) of the
        #    statement i == mutant_aa
        aa_bool.append( i == mutant_aa )

    # modify the mutating residue's assignment in the PackerTask using the
    #    Vector1 of booleans across the proteogenic amino acids
    task.nonconst_residue_task( mutant_position
        ).restrict_absent_canonical_aas( aa_bool )

    # prevent residues from packing by setting the per-residue "options" of
    #    the PackerTask
    center = pose.residue( mutant_position ).nbr_atom_xyz()
    for i in range( 1 , pose.total_residue() + 1 ):
        # only pack the mutating residue and any within the pack_radius
        if not i == mutant_position or center.distance_squared(
                test_pose.residue( i ).nbr_atom_xyz() ) > pack_radius**2:
            task.nonconst_residue_task( i ).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover( pack_scorefxn , task )
    packer.apply( test_pose )

    return test_pose

###############################################################################

## Protonate selected residues in the pose (for Moon & Fleming case only)
def adjust_for_pH( pose, resnum ):

    ## Create a copy of the pose to return
    adjusted_pose = Pose()
    adjusted_pose.assign( pose )

    ## Apply residue variant type at the mutated position
    add_variant_type_to_pose_residue( adjusted_pose, resnum, PROTONATED )

    # Create a membrane framework pH scorefunction
    pH_scorefxn = create_score_function( "fa_menv_pHmode_2014" )

    # Create a packer task (standard)
    task = TaskFactory.create_packer_task( adjusted_pose )

    # Prevent all residues from repacking except for the 
    # protonated residue by setting per-residue options
    # of the packer task
    for i in range( 1, adjusted_pose.total_residue() + 1 ): 
        task.nonconst_residue_task( i ).prevent_repacking()

    # Repack the protonated residue, restricting to protonated rotamers
    # only
    packer = PackRotamersMover( pH_scorefxn , task )
    packer.apply( adjusted_pose )

    return adjusted_pose


if __name__ == "__main__" : main(sys.argv)

