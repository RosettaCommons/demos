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

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/types.hh>
#include <core/scoring/sasa.hh>
#include <core/util/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/keys/OptionKeys.hh>

#include <core/io/database/open.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

using namespace ObjexxFCL::fmt;

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] )
{
	// options, random initialization
	core::init( argc, argv );

	using namespace core::options;
	using namespace core::options::OptionKeys;

	// std::string pdbfile = option[ phil::s ]();
	utility::vector1< std::string > pdbfiles = option[ james::pdbfile ]();

	utility::vector1< std::string >::const_iterator iter;
	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter ) {
		std::string pdbfile = *iter;
		if ( pdbfile == "" ) {
			utility_exit_with_message( "Unable to open file: " + pdbfile + '\n' );
		}
		core::pose::PoseOP mypose ( new core::pose::Pose );
		std::cerr << "READING " << pdbfile << std::endl;
		core::io::pdb::pose_from_pdb( *mypose, pdbfile ); // default is standard fullatom residue_set

		std::string outfile  = pdbfile + ".distance_map";
		std::ofstream output( outfile.c_str() );
		if ( ! output.is_open() ) {
			utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
		}

		// output << A( 10, "index ");
		// for ( unsigned int i = 1; i <= mypose->total_residue(); ++i ) {
		// 	output << I(10, i);
		// }
		// output << std::endl;

		for ( unsigned int i = 1; i <= mypose->total_residue(); ++i ) {
			// output << I( 10, i );
			for ( unsigned int j = 1; j <= mypose->total_residue(); ++j ) {
				core::conformation::Residue resi = mypose->residue(i);
				core::conformation::Residue resj = mypose->residue(j);

				int atomno = 2; // c-alpha atoms
				core::Real distance  = mypose->residue(i).xyz(atomno).distance( resj.xyz(atomno) );
				output 	<< F( 10, 4, distance );

			} // 	for ( unsigned int j = i + 1; j <= mypose->total_residue(); ++j )
			output << std::endl;
		}		// for ( unsigned int i = 1; i <= mypose->total_residue(); ++i )
		output.close();
	} // 	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter )
} // int main( int argc, char * argv [] )
