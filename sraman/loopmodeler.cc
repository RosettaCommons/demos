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
#include <devel/IntegratedLoop/LoopModeler.hh>
#include <protocols/loops/LoopClass.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

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
using namespace fragment;


void* my_main( void * )
{

	using namespace core::options;
	using namespace core::options::OptionKeys;

	std::string native_fn   = option[ in::file::native ]();
	std::string frag3_file  = option[ in::file::frag3 ]();
	std::string frag9_file  = option[ in::file::frag9 ]();
	std::string loop_file = option[LoopModel::loop_file]();


	core::pose::Pose native_pose;
	io::pdb::pose_from_pdb( native_pose, native_fn );

	std::string sequence = native_pose.sequence(); // must match sequence of fragments! maybe read from fasta instead

	ConstantLengthFragSetOP fragset3mer = new ConstantLengthFragSet( 3 );
	ConstantLengthFragSetOP fragset9mer = new ConstantLengthFragSet( 9 );
	fragset3mer->read_fragment_file( frag3_file );
	fragset9mer->read_fragment_file( frag9_file );

	//Loops file
	protocols::loops::Loops loop_list;
	loop_list.read_file( loop_file );

	// pose::Pose
	pose::Pose input_pose;
	//core::io::pdb::pose_from_pdb this will come here

	std::cout << "LoopModeler start" << std::endl;

	// create list of tags already processed

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

	for( Size iteration = 1; iteration <= (Size) nstruct_flag; iteration++ ){
		std::string output_tag = "S_" + right_string_of(iteration, 8, '0');
		if ( tags_done[ output_tag ] ) {
			std::cerr << "Tag: " << output_tag << " - already processed" << std::endl;
			continue;
		}

		std::cerr << "Tag: " << output_tag << " - processing." << std::endl;

		pose::PoseOP loop_model_pose = new pose::Pose( input_pose );


		devel::LoopModel::LoopModeler myLoopModeler( loop_list, fragset3mer, fragset9mer );
		//		myLoopModeler.init( *loop_model_pose ); // Do I need this ?


		core::scoring::ScoreFunctionOP scorefxn
			= core::scoring::ScoreFunctionFactory::create_score_function( "score3" );

		int start_time = time(NULL);
		//		myLoopModeler.apply( *fold_pose ); //Uncomment this when everything is done
		int end_time   = time(NULL);
		std::cout << "TIMEFORLOOPMODELING: " << end_time - start_time << std::endl;

		bool fullatom( false );

		std::string silent_file = option[ out::file::silent ]();
		if ( option[ out::file::silent ].user() ) {
			core::io::silent::ProteinSilentFileData sfd;
			(*scorefxn)(*loop_model_pose);
			(*scorefxn).show( std::cout, *loop_model_pose );
			sfd.write_pose( *loop_model_pose, native_pose, scorefxn, silent_file, output_tag, fullatom );
		} else {
			std::cout << "Warning: You are evil for outputting PDBs! Use -silent_out instead.\n";
			std::string filename =  output_tag + ".pdb";
			io::pdb::dump_pdb( *loop_model_pose , filename );
		}
	}
}



int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );
	protocols::viewer::viewer_main( my_main );
	return 0;
}

