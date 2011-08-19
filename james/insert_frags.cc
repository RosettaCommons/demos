// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers


#include <core/types.hh>
#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/Frame.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <utility/options/FileVectorOption.hh>

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

	// pose.conformation().insert_chain_ending( pose.total_residue() - 1 );		// probably not necessary

} // make_pose_match_sequence_

void* my_main( void * )
{

	using namespace core::options;
	using namespace core::options::OptionKeys;

	std::string native_fn   = option[ in::file::native ]();
	utility::vector1< int > frag_sizes( option[ in::file::frag_sizes ]() );
	FileVectorOption frag_files( option[ in::file::frag_files ] );

	core::pose::Pose native_pose;
	io::pdb::pose_from_pdb( native_pose, native_fn );

	Size frag_length = frag_sizes[1];
	ConstantLengthFragSetOP fragset = new ConstantLengthFragSet( frag_length );
	fragset->read_fragment_file( frag_files[1] );

	pose::Pose current_pose, old_pose;
	std::string sequence = native_pose.sequence();
	make_pose_from_sequence_( sequence,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )),
		current_pose
	);

	// make extended chain
	for ( Size pos = 1; pos <= current_pose.total_residue(); pos++ ) {
		current_pose.set_phi( pos, -45 );
		current_pose.set_psi( pos, -45 );
		current_pose.set_omega( pos, 180 );
	}

	// make a MoveMap
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	for ( Size i = 1; i <= current_pose.total_residue(); ++i ) {
		movemap->set_bb( true );
	}

	// make a scoring function
	core::scoring::ScoreFunctionOP scorefxn(
		new core::scoring::ScoreFunction()
	);
	scorefxn->set_weight( core::scoring::vdw, 1.0 );

	if ( option[ james::debug ]() ) {
		current_pose = native_pose;
	}
	core::io::silent::ProteinSilentFileData sfd;
	old_pose = current_pose; // used to roll back moves!

	/// BEGIN SINGLE-FRAGMENT INSERTIONS!!!
	std::string outfile = "single." + option[ out::file::silent ]();
	core::fragment::Frame frame;
	Size frag_size = fragset->max_frag_length();
	Size decoy_count = 1;

	Real vdw_start = (*scorefxn)(current_pose);
	std::cout << "vdw_start = " << vdw_start << std::endl;

	std::cout << "frag_size = " << frag_size << std::endl;
	for ( Size pos = 1; pos < current_pose.total_residue() - frag_size + 1; ++pos ) {
		bool check1 = fragset->get_frame( pos, frame );
		if ( !check1 ) {
			std::cerr << "Warning: can't find frame at position " << pos << std::endl;
			continue;
		}

		for ( Size nfrag = 1; nfrag <= frame.nr_frags(); ++nfrag ) {
			frame.apply( *movemap, nfrag, current_pose );

			Real vdw_delta = vdw_start - (*scorefxn)(current_pose);
			core::Real CA_rmsd = core::scoring::CA_rmsd( old_pose, current_pose );

			std::string decoy_tag = "S_" + string_of( decoy_count, 12 );

			core::io::silent::ProteinSilentStruct silent_struct(
				current_pose, decoy_tag, false
			);

			silent_struct.clear_energies();
			silent_struct.add_energy( "CA_rmsd",   CA_rmsd    );
			silent_struct.add_energy( "vdw_delta", vdw_delta  );
			silent_struct.add_energy( "pos",       pos        );
			silent_struct.add_energy( "nfrag",     nfrag      );

			sfd.write_silent_struct( silent_struct, outfile  );
			++decoy_count;
		}

		current_pose = old_pose;
	}

	// BEGIN DOUBLE-FRAGMENT INSERTIONS!!!
	outfile = "double." + option[ out::file::silent ]();

	core::fragment::Frame frame1, frame2;
	decoy_count = 1;
	for ( Size pos1 = 1; pos1 < current_pose.total_residue() - 2 * frag_size + 1; ++pos1 ) {
		bool check1 = fragset->get_frame( pos1, frame1 );
		if ( !check1 ) {
			std::cerr << "Warning: can't find frame at position " << pos1 << std::endl;
			continue;
		}

		Size pos2 = pos1 + frag_size;
		std::cout << "inserting frags at (" << pos1 << "," << pos2 << ")" << std::endl;

		bool check2 = fragset->get_frame( pos2, frame2 );
		if ( !check2 ) {
			std::cerr << "Warning: can't find frame at position " << pos2 << std::endl;
			continue;
		}

		for ( Size nfrag1 = 1; nfrag1 <= frame1.nr_frags(); ++nfrag1 ) {
			frame1.apply( *movemap, nfrag1, current_pose );

			for ( Size nfrag2 = 1; nfrag2 <= frame2.nr_frags(); ++nfrag2 ) {
				frame2.apply( *movemap, nfrag2, current_pose );

				Real vdw_delta = vdw_start - (*scorefxn)(current_pose);
				core::Real CA_rmsd = core::scoring::CA_rmsd( old_pose, current_pose );

				std::string decoy_tag = "S_" + string_of( decoy_count, 12 );
				core::io::silent::ProteinSilentStruct silent_struct(
					current_pose, decoy_tag, false
				);

				silent_struct.clear_energies();
				silent_struct.add_energy( "CA_rmsd",   CA_rmsd   );
				silent_struct.add_energy( "vdw_delta", vdw_delta );
				silent_struct.add_energy( "pos1",      pos1      );
				silent_struct.add_energy( "pos2",      pos2      );
				silent_struct.add_energy( "nfrag1",    nfrag1    );
				silent_struct.add_energy( "nfrag2",    nfrag2    );

				sfd.write_silent_struct( silent_struct, outfile  );
				++decoy_count;
			} // nfrag1
		} // nfrag2

		current_pose = old_pose;
	} // for ( it1 = fragset->begin(), it_end = fragset->end(); it1 != it_end; ++it1; )

	return NULL;
}

int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );
	protocols::viewer::viewer_main( my_main );
	return 0;
}
