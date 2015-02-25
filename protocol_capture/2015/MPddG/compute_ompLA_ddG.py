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
import sys, os
import commands

# rosetta
import rosetta.protocols.membrane

## Default Values
repack_radius = 0.0


import random
from rosetta import Pose
from rosetta import create_score_function
from rosetta import TaskFactory
from rosetta.utility import vector1_bool
from rosetta import aa_from_oneletter_code
from rosetta import PackRotamersMover
from rosetta.core.pose import PDBInfo

## @brief Compute ddG of mutation in a protein at specified residue and AA position
def compute_ddG( pose, sfxn, resnum, aa ): 

    # Score Native Pose
    native_score = sfxn( pose )

    # Perform Mutation at residue <resnum> to amino accid <aa>
    mutated_pose = mutate_residue( pose, resnum, aa, repack_radius, sfxn )

    # Adjust for pH
    #adjust_for_pH( pose, resnum )

    # Score Mutated Pose
    mutant_score = sfxn( mutated_pose )

    # Print resulting ddG
    return mutant_score - native_score

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

## @brief Main - Add Membrane to Pose, Compute ddG
def main( argv ):

    rosetta.init( extra_options="-mp:setup:spanfiles inputs/1qd6_tr_C.span -run:constant_seed -in:ignore_unrecognized_res" )

    # Load Pose, & turn on the membrane
    pose = pose_from_pdb( "inputs/1qd6_tr_C.pdb" );
    sfxn = create_score_function( "mpframework_smooth_fa_2014" );

    #sfxn = create_score_function( "mpframework_pHmode_fa_2014" );
    
    # Add Membrane to Pose
    add_memb = rosetta.protocols.membrane.AddMembraneMover()
    add_memb.apply( pose )

    # Setup in a topology based membrane
    init_mem_pos = rosetta.protocols.membrane.MembranePositionFromTopologyMover()
    init_mem_pos.apply( pose )

    # Compute mutations 
    #AAs = ['D', 'E']
    AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    ddGs = []
    for aa in AAs:
        print compute_ddG( pose, sfxn, 181, aa  )

if __name__ == "__main__" : main(sys.argv)

