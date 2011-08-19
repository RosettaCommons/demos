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
/// @author James Thompson
/// @author Greg Taylor

// libRosetta headers

#include <core/options/keys/OptionKeys.hh>
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/sequtil.hh>

#include <core/util/basic.hh>
#include <core/util/Tracer.hh>
#include <core/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#include <utility/file/FileName.hh>

#include <protocols/viewer/viewers.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/util/Tracer.hh>

// GKT - added
#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>


using core::util::T;
using core::util::Warning;
using core::util::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

using namespace core;
using namespace fragment; // GKT
using namespace protocols; // GKT
using namespace abinitio; // GKT

#include "homolog_cst.hh"

///////////////////////////////////////////////////////////////////////////////
void*
build_template( void* )
{
	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP rsd_set
		= core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	core::pose::Pose native_pose, fold_pose;
	if ( option[ in::file::native ].user() ) {
		std::string native_fn = option[ in::file::native  ]();
		io::pdb::pose_from_pdb( native_pose, native_fn );
	}

	std::string sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]() )[1];

	chemical::make_pose_from_sequence( fold_pose, sequence, *rsd_set );
	add_constraints( &fold_pose, option[ OptionKeys::constraints::cst_file ]() );

	fold_pose.dump_pdb( "starting_template.pdb" );

	using namespace optimization;
	kinematics::MoveMap mm;

	core::scoring::ScoreFunctionOP bump_scorefxn( new core::scoring::ScoreFunction() );
	bump_scorefxn->set_weight( vdw, 1.0 );
	bump_scorefxn->set_weight( atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ]() );

	mm.set_bb ( true );

	fold_pose.dump_pdb("before_fold.pdb");

	//AtomTreeMinimizer().run(
	//		fold_pose, mm, (*bump_scorefxn), MinimizerOptions( "dfpmin_armijo_nonmonotone", 0.001, true )
	//	);
	//fold_pose.dump_pdb("denatured_fold.pdb");

	// Read the fragment file for 3 and 9 aa fragments
	ConstantLengthFragSetOP fragset3mer = new ConstantLengthFragSet( 3 );
	ConstantLengthFragSetOP fragset9mer = new ConstantLengthFragSet( 9 );
	fragset3mer->read_fragment_file( option[ in::file::frag3 ]() );
	fragset9mer->read_fragment_file( option[ in::file::frag9 ]() );

	// make a MoveMap - define during an optimization which aa's are allowed to move
	// for this protocol all amino acids should be allowed to freely move.
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	for ( Size i = 1; i <= fold_pose.total_residue(); ++i ) {
		movemap->set_bb( true );
	}

	FoldConstraints abinitio_protocol( fragset3mer, fragset9mer, movemap );
	abinitio_protocol.init( fold_pose );
	abinitio_protocol.set_constraint_weight( option[ OptionKeys::constraints::cst_weight]() );

	scoring::ScoreFunctionOP scorefxn
			= scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scorefxn->set_weight( scoring::atom_pair_constraint, option[ OptionKeys::constraints::cst_weight ]() );

	abinitio_protocol.apply( fold_pose );

	// write the silent files
	core::io::silent::ProteinSilentFileData sfd;
	sfd.write_pose( fold_pose, "silent_file", "whatever", false );

	fold_pose.dump_pdb("after_fold.pdb");

	return 0;
} // build_template

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );
	protocols::viewer::viewer_main( build_template );

	return 0;
}
