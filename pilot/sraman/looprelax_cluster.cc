// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file demo/sraman/looprelax.cc
/// @brief demo program for implementing looprelax protocol
/// @author Srivatsan Raman
/// @author James Thompson


// mini headers
#include <core/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/util/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/options/util.hh>
#include <core/options/after_opts.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>
#include <numeric/random/random.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <devel/loops/looprelax_protocols.hh>
#include <protocols/viewer/viewers.hh>

// C++ headers
#include <iostream>
#include <string>
/*
int
main( int argc, char * argv [] )
{

	core::init(argc, argv);

	//	numeric::random::RandomGenerator::initializeRandomGenerators(
	//		 111111, numeric::random::_RND_TestRun_, "ran3");

	using namespace protocols::moves;
	using namespace scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	//	std::string pdbfile = start_file(); //option[ phil::s ]();
	std::string pdbfile = option[ OptionKeys::looprelax::input_pdb ]().name();
	pose::PoseOP pose ( new pose::Pose );
	io::pdb::pose_from_pdb( *pose, pdbfile ); // defaults to fa_standard residues

	//	std::string wts_tag = STANDARD_WTS;
	//	core::scoring::ScoreFunctionOP scorefxn(
	//		ScoreFunctionFactory::create_score_function( wts_tag )

	devel::looprelax::LoopRelax::loop_relax_main( pose );
//	looprelax.apply();
}
*/
void
loop_relax_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::scoring;

	core::pose::Pose pose;
	std::cout << "loop_relax_viewer" << std::endl;
	core::io::pdb::pose_from_pdb( pose, option[ OptionKeys::looprelax::input_pdb ]().name() );
	std::string silent_filename = option[ out::file::silent ]();

	core::io::silent::ProteinSilentFileData sfd( silent_filename );
	core::io::silent::Structure_Map::const_iterator iter;
	utility::vector1< std::string > tags = sfd.read_tags_fast( silent_filename );
	utility::vector1< std::string >::const_iterator i;
	std::map< std::string, bool > tags_done;
	for ( i = tags.begin(); i != tags.end(); ++i ) {
		tags_done[ *i ] = true;
	}

	int const nstruct_flag = option[ out::nstruct ];

	for( int iteration = 1; iteration <= nstruct_flag; iteration++ ){
		std::string output_tag = "S_" + right_string_of(iteration, 8, '0');
		if ( tags_done[ output_tag ] ) {
			std::cerr << "Tag: " << output_tag << " - already processed" << std::endl;
			continue;
		}
		std::cerr << "Tag: " << output_tag << " - processing." << std::endl;

	devel::looprelax::LoopRelax myLoopRelax;

	int start_time = time(NULL);
	myLoopRelax.apply( pose );
	int end_time = time( NULL );
	std::cout << "TIME FOR LOOPRELAX " << end_time - start_time << std::endl;

	//scoring
	std::string function_tag("cen_std"),patch_tag("score4L");
	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( function_tag, patch_tag ) );

	bool fullatom( false );
	std::string silent_file = option[ out::file::silent ]();
	core::io::silent::ProteinSilentFileData sfd = new core::io::silent::ProteinSilentFileData();
	if ( option[ out::file::silent ].user() ) {
		(*scorefxn)( pose );
		sfd.write_pose( pose, silent_file, output_tag, fullatom );
	} else {
		std::cout << "Warning: You are evil for outputting PDBs! Use -silent_out instead.\n";
		std::string filename =  output_tag + ".pdb";
		core::io::pdb::dump_pdb( pose , filename );
	}
	}
}



void*
my_main( void* )
{
	numeric::random::RandomGenerator::initializeRandomGenerators(
		 111111, numeric::random::_RND_TestRun_, "ran3");

	loop_relax_test();
	return 0;
}

int
main( int argc, char * argv [] )
{

	// initialize option and random number system
	core::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}
