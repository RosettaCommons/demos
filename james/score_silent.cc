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

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/keys/OptionKeys.hh>

#include <core/util/basic.hh>
#include <core/util/Tracer.hh>
#include <core/io/database/open.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>

#include <utility/vector1.hh>

#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Warning;
using core::util::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

int
main( int argc, char* argv [] )
{
	// options, random initialization
	core::init( argc, argv );

	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	// configure score function
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();

	// configure silent-file data object
	core::io::silent::ProteinSilentFileData sfd;

	std::string infile  = option[ in ::file::silent ]();
	std::string outfile = option[ out::file::silent ]();

	if ( option[ in::file::silent ].user() ) {
		sfd.read_file( infile );
	}

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		core::pose::Pose pose;
		iter->fill_pose( pose, *rsd_set );
		(*scorefxn)(pose);
		sfd.write_pose( pose, outfile, iter->decoy_tag(), option[ in::file::fullatom ]() );
	}

	return 0;
}
