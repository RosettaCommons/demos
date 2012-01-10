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

#include <numeric/random/random.hh>

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

	core::pose::Pose pose;
	std::cout << "loop_relax_viewer" << std::endl;
	core::io::pdb::pose_from_pdb( pose, option[ OptionKeys::looprelax::input_pdb ]().name() );
	devel::looprelax::LoopRelax myLoopRelax;
	myLoopRelax.apply( pose );

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
