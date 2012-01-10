// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file demo/rms_demo.cc
/// @brief

/// @author James Thompson

#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>
#include <core/types.hh>

#include <core/init.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pose_io.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );
	using namespace core::options;
	using namespace core::options::OptionKeys;

	utility::vector1< std::string > pdbfiles = option[ james::pdbfile ]();

	std::string file1 = pdbfiles[1];
	std::string file2 = pdbfiles[2];

	core::pose::PoseOP pose1 ( new core::pose::Pose );
	core::pose::PoseOP pose2 ( new core::pose::Pose );

	core::io::pdb::pose_from_pdb( *pose1, file1 );
	core::io::pdb::pose_from_pdb( *pose2, file2 );

	std::cout << "Calculating RMSD ... " << std::endl;
	if ( pose1->total_residue() == pose2->total_residue() ) {
		std::cout << "all-atom RMSD = " << core::scoring::all_atom_rmsd( *pose1, *pose2 ) << std::endl;
		std::cout << "CA_rmsd = " << core::scoring::CA_rmsd( *pose1, *pose2 ) << std::endl;
		std::cout << "CA_maxsub = " << core::scoring::CA_maxsub( *pose1, *pose2 ) << std::endl;
	} else {
		std::cerr
			<< "Error: can't currently calculate RMSDs between structures with different lengths."
			<< std::endl;
	}

	return 0;
}
