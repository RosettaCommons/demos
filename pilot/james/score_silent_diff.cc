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
#include <core/scoring/ScoringManager.hh>
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
// #include <core/io/silent/EmptySilentStruct.hh>
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
	} else if ( option[ james::pdbfile ].user() ) {
		core::pose::Pose temp_pose;
		utility::vector1< std::string > pdbfiles = option[ james::pdbfile ]();
		for ( utility::vector1< std::string >::iterator it = pdbfiles.begin(), end = pdbfiles.end(); it != end; ++it ) {
			core::io::pdb::pose_from_pdb( temp_pose, *rsd_set, *it );
			core::io::silent::ProteinSilentStructOP ss(
				new core::io::silent::ProteinSilentStruct( temp_pose, *it, option[ in::file::fullatom ]() )
			);
			sfd.add_structure( ss );
		}
	} else {
		utility_exit_with_message( "Error: can't get any structures! Use -in::file::silent or -james::pdbfile!" );
	}

	utility::vector1< std::string > energies_wanted;
	energies_wanted.push_back( "env"     );
	energies_wanted.push_back( "pair"    );
	energies_wanted.push_back( "cbeta"   );
	energies_wanted.push_back( "rg"      );
	energies_wanted.push_back( "vdw"     );
	energies_wanted.push_back( "hs_pair" );
	energies_wanted.push_back( "ss_pair" );
	energies_wanted.push_back( "rsigma"  );

	std::map< std::string, utility::vector1< core::Real > > energy_diffs;

	core::pose::Pose pose;
	sfd.set_fullatom( option[ in::file::fullatom ]() );
	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		iter->fill_pose( pose, *rsd_set );

		(*scorefxn)(pose);
		scorefxn->accumulate_residue_total_energies( pose );
		core::io::silent::ProteinSilentStruct ss( pose, iter->decoy_tag(), option[ in::file::fullatom ]() );

		std::map< std::string, core::Real > energies;
		for ( utility::vector1< std::string >::const_iterator it = energies_wanted.begin(), end = energies_wanted.end();
					it != end; ++it ) {

			core::Real old_e  = iter->get_energy( *it );
			core::Real new_e  =    ss.get_energy( *it );

			// ss.add_energy( *it + "_old", old_e );
			// ss.add_energy( *it + "_new", new_e );
			// ss.add_energy( *it + "_diff", new_e - old_e );
			// energies[ *it + "_old"  ] = old_e;
			// energies[ *it + "_new"  ] = new_e;
			energies[ *it + "_diff" ] = new_e - old_e;
		}
		ss.clear_energies();
		for ( std::map< std::string, core::Real >::const_iterator it = energies.begin(), end = energies.end();
					it != end; ++it ) {
			ss.add_energy( it->first, it->second );
		}

		sfd.write_silent_struct( ss, outfile );
	}

	return 0;
}
