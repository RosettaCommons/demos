// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// libRosetta headers


#include <core/types.hh>
#include <core/init.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/viewer/viewers.hh>
#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;



using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;

void
make_pose_from_sequence_(
	std::string sequence,
	chemical::ResidueTypeSet const& residue_set,
	pose::Pose& pose
) {
	using namespace chemical;
	// clear all of the old data in the pose
	pose.clear();

	// setup the pose by appending the appropriate residues residues
	for ( Size seqpos = 1; seqpos <= sequence.length(); ++seqpos ) {
		char aa = sequence[seqpos-1]; // string indexing is zero-based!
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueTypeCAPs const & rsd_type_list( residue_set.aa_map( my_aa ) );
		Size best_index = 1;
		ResidueType const & rsd_type( *(rsd_type_list[ best_index ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		if ( seqpos == 1 ) {
			pose.append_residue_by_jump( *new_rsd, 1 );
		} else {
			pose.append_residue_by_bond( *new_rsd, true );
		}
	} // for seqpos

	pose.conformation().insert_chain_ending( pose.total_residue() - 1 );		// probably not necessary

} // make_pose_match_sequence_

//int main (int argc, char* argv[])
void* my_main( void * )
{

	//  core::init(argc, argv);

	std::string sequence( "MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE" ); //GB3


	ConstantLengthFragSetOP fragset3mer = new ConstantLengthFragSet( 3 );
	fragset3mer->read_fragment_file( "input_abinitio/mfr_aa2GB3_03_05.200_v1_3" );
	ConstantLengthFragSetOP fragset9mer = new ConstantLengthFragSet( 9 );
	fragset9mer->read_fragment_file( "input_abinitio/mfr_aa2GB3_09_05.200_v1_3" );


	pose::Pose extended_pose;
	make_pose_from_sequence_( sequence,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )),
		extended_pose
	);

	std::string const pdbfile ( "input_abinitio/2GB3.pdb" );
	std::cout << "coarse start" << std::endl;
	// io::pdb::pose_from_pdb( extended_pose, pdbfile );

	// make extended chain
	for ( Size pos = 1; pos<= extended_pose.total_residue(); pos++ ) {
		extended_pose.set_phi( pos, -45 );
		extended_pose.set_psi( pos, -45 );
		extended_pose.set_omega( pos, 180 );
	}

	pose::Pose fold_pose = extended_pose;
	io::pdb::dump_pdb( extended_pose , "starting_pose.pdb" );

	kinematics::MoveMapOP movemap = new kinematics::MoveMap; // not used yet
	ClassicAbinitio abinitio( fragset3mer, fragset9mer, movemap );
	abinitio.init( fold_pose );
	// protocols::viewer::add_monte_carlo_viewer( abinitio.mc() );

	abinitio.apply( fold_pose );
	io::pdb::dump_pdb( fold_pose , "resulting_pose.pdb" );

	return NULL;
}


int
main( int argc, char * argv [] )
{
	// options, random initialization
	std::cerr<< "Initialize rosetta core... " << std::endl;
	core::init( argc, argv );
	std::cerr<< "GLUT call viewer_main " << std::endl;
	protocols::viewer::viewer_main( my_main );
	return 0;
}

