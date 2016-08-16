#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file: predict_ompLA_ddG.py
##
## @brief: Predict the ddG of mutation at position 181 of OmpLA
## @details: Use the rosetta membrane framework (RosettaMP) to predict the ddG of 
## mutation from alanine to all canonical residues at poisiton 181 of OmpLA. This
## scripts uses the mutate_residue function, adapted from the mutate.py script in the 
## 2012 PyRosetta workshops written by Evan Baugh. 
##
## Citation: Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from 
## transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. Available: 
## http://www.pnas.org/content/early/2011/05/16/1103979108.
##
## @author: Rebecca F. Alford (rfalford12@gmail.com)

## *** Note - this script is not currently using the pH correction ** ##

from rosetta import *

# tools
import sys, os
import commands
import random

# rosetta
import rosetta.protocols.membrane
from rosetta import Pose
from rosetta import create_score_function
from rosetta import TaskFactory
from rosetta.utility import vector1_bool
from rosetta import aa_from_oneletter_code
from rosetta import PackRotamersMover
from rosetta.core.pose import PDBInfo

###############################################################################
## @brief:  Main Function for ddG Prediction
## @details: Load in pose, setup membrane energy function, setup the membrane
## framework, calculate ddGs 
def main( argv ):

    ## Step 1: Initialize Rosetta, including spanfile for the membrane framework
    rosetta.init( extra_options="-mp:setup:spanfiles inputs/1qd6_tr_C.span -run:constant_seed" )

    ## Step 2: Load in pose as a membrane pose (using AddMembraneMover)
    pose = pose_from_pdb( "inputs/1qd6_tr_C.pdb" );
    add_memb = rosetta.protocols.membrane.AddMembraneMover()
    add_memb.apply( pose )
    
    ## Step 3: Orient Pose in the membrane based on transmembrane spans 
    init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
    init_mem_pos.apply( pose )

    ## Step 4: Create a membrane energy function
    ## @note: Two possible energy functions: (1) standard membrane full atom energy function
    ## or (2) energy function accommodating for pH (requires extra flags - see README)
    sfxn = create_score_function( "mpframework_smooth_fa_2012" );

    ## Step 5: Repack neighbors within 8.0 angstroms of the mutant position
    repacked_native = mutate_residue( pose, 181, 'A', 8.0, sfxn )

    ## Step 6: Compute the ddG of mutation from alanine to each of the 20 canonical
    ## amino acids. Then print the ddG file to an output file (ompLA_ddG.out)
    AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    with file( 'ompLA_ddG.out', 'a' ) as f: 
        f.write( "Res AA ddG\n" )
        for aa in AAs: 
            ddG = compute_ddG( repacked_native, sfxn, 181, aa, 8.0 )
            f.write( str(181) + " " + aa + " " + str( round( ddG, 3 ) ) + "\n" )
    f.close

###############################################################################
## @brief Compute the ddG of mutation
## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, sfxn, resnum, aa, repack_radius ): 

    # Score Native Pose
    native_score = sfxn( pose )

    # Perform Mutation at residue <resnum> to amino accid <aa>
    mutated_pose = mutate_residue( pose, resnum, aa, repack_radius, sfxn )

    # Score Mutated Pose
    mutant_score = sfxn( mutated_pose )

    # Print resulting ddG
    return mutant_score - native_score

###############################################################################
# @brief Replace the residue at <resid> in <pose> with <new_res> and allows
# repacking within a given <pack_radius> 
def mutate_residue( pose, mutant_position, mutant_aa, pack_radius, sfxn ):

    if pose.is_fullatom() == False:
        IOError( 'mutate_residue only works with fullatom poses' )

    test_pose = Pose()
    test_pose.assign( pose )

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
    for i in range( 1, pose.total_residue() + 1 ): 
        dist = center.distance_squared( test_pose.residue( i ).nbr_atom_xyz() );  
        # only pack the mutating residue and any within the pack_radius
        if i != mutant_position and dist > pow( float( pack_radius ), 2 ) :
            task.nonconst_residue_task( i ).prevent_repacking()

    # apply the mutation and pack nearby residues
    packer = PackRotamersMover( sfxn , task )
    packer.apply( test_pose )

    return test_pose

if __name__ == "__main__" : main(sys.argv)

