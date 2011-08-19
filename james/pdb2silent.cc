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

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <protocols/loops/ccd_closure.hh>
#include <protocols/relax_protocols.hh>

#include <core/options/keys/OptionKeys.hh>

#include <core/types.hh>

#include <core/util/prof.hh> // profiling

#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/util/basic.hh>
#include <core/util/Tracer.hh>
#include <core/io/database/open.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/random/ran3.hh>

#include "homolog_cst.hh"

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

	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	utility::vector1< std::string > pdbfiles = option[ james::pdbfile ]();
	std::string silent_outfile = option[ out::file::silent ]();
	core::pose::Pose original_pose, rebuilt_pose;
	core::io::silent::ProteinSilentFileData sfd;

	sfd.set_fullatom( option[ in::file::fullatom ]() );
	for ( utility::vector1< std::string >::iterator it = pdbfiles.begin(), end = pdbfiles.end(); it != end; ++it ) {
		core::io::pdb::pose_from_pdb( original_pose, *rsd_set, *it );

		utility::file::file_delete( silent_outfile );
		sfd.write_pose( original_pose, silent_outfile, *it + ".original", option[ in::file::fullatom ]() );
		core::io::silent::ProteinSilentStruct pss( original_pose, "original", option[ in::file::fullatom ]() );
		pss.print_conformation( std::cout );

		for ( Size i = 1; i <= original_pose.total_residue(); ++i ) {
			std::cout << I( 8, i )
								<< F( 8, 3, original_pose.phi(i)   )
								<< F( 8, 3, original_pose.psi(i)   )
								<< F( 8, 3, original_pose.omega(i) )
								<< "\n";
		}

		for ( int ntries = 1; ntries <= 10; ++ntries ) {
			sfd.read_file( silent_outfile );
			// utility::file::file_delete( silent_outfile );

			core::io::silent::ProteinSilentFileData::iterator iter = sfd.begin();
			for ( iter = sfd.begin(); iter != sfd.end(); ++iter ) {
				iter->fill_pose( rebuilt_pose, *rsd_set );
				iter->print_conformation( std::cout );
			}
			if ( iter == sfd.end() ) {
				std::cout << "True!!!" << std::endl;
			} else {
				std::cout << "False!!!" << std::endl;
			}

			std::cout << "rmsd = " << core::scoring::CA_rmsd( original_pose, rebuilt_pose ) << "\n";
			sfd.write_pose( rebuilt_pose, silent_outfile, *it + "." + string_of(ntries), option[ in::file::fullatom ]() );
		}
	}

	return 0;
}
