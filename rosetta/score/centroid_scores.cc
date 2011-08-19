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
//#include <core/options/option.hh>

#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/VariantType.hh>

#include <core/init.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>



#include "utility/pointer/owning_ptr.hh"
#include <utility/vector1.hh>

#include "ObjexxFCL/formatted.o.hh"


// C++ headers
#include <iostream>
#include <string>
#include <vector>


void
print_scores( )
{

	using namespace core;
	using namespace chemical;
	using namespace pose;
	using namespace ObjexxFCL::fmt;
	using namespace scoring;
	using namespace utility;
	using std::vector;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::string;
	using core::Real;

	vector1<string> pdbs;
	pdbs.push_back("1a19"); pdbs.push_back("1a32"); pdbs.push_back("1a68"); pdbs.push_back("1acf");
	pdbs.push_back("1ail"); pdbs.push_back("1aiu"); pdbs.push_back("1b3a"); pdbs.push_back("1bgf");
	pdbs.push_back("1bk2"); pdbs.push_back("1bkr"); pdbs.push_back("1bm8"); pdbs.push_back("1bq9");
	pdbs.push_back("1c8c"); pdbs.push_back("1c9o"); pdbs.push_back("1cc8"); pdbs.push_back("1cei");
	pdbs.push_back("1cg5"); pdbs.push_back("1ctf"); pdbs.push_back("1dhn"); pdbs.push_back("1e6i");
	pdbs.push_back("1elw"); pdbs.push_back("1enh"); pdbs.push_back("1ew4"); pdbs.push_back("1eyv");
	pdbs.push_back("1fkb"); pdbs.push_back("1fna"); pdbs.push_back("1gvp"); pdbs.push_back("1hz6");
	pdbs.push_back("1ig5"); pdbs.push_back("1iib"); pdbs.push_back("1kpe"); pdbs.push_back("1lis");
	pdbs.push_back("1lou"); pdbs.push_back("1nps"); pdbs.push_back("1opd"); pdbs.push_back("1pgx");
	pdbs.push_back("1ptq"); pdbs.push_back("1r69"); pdbs.push_back("1rnb"); pdbs.push_back("1scj");
	pdbs.push_back("1shf"); pdbs.push_back("1ten"); pdbs.push_back("1tif"); pdbs.push_back("1tig");
	pdbs.push_back("1tit"); pdbs.push_back("1tul"); pdbs.push_back("1ubi"); pdbs.push_back("1ugh");
	pdbs.push_back("1urn"); pdbs.push_back("1utg"); pdbs.push_back("1vcc"); pdbs.push_back("1vie");
	pdbs.push_back("1vls"); pdbs.push_back("1who"); pdbs.push_back("1wit"); pdbs.push_back("256b");
	pdbs.push_back("2acy"); pdbs.push_back("2chf"); pdbs.push_back("2ci2"); pdbs.push_back("2vik");
	pdbs.push_back("4ubp"); pdbs.push_back("5cro");

	std::cout << "SCORE "
						<< " pdb "
						<< "       env"
						<< "      pair"
						<< "       vdw"
						<< "        hs"
						<< "        ss"
						<< "     sheet"
						<< "     cbeta"
						<< "    rsigma"
						<< "   hb_srbb"
						<< "   hb_lrbb"
						<< "        rg"
						<< "      rama"
						<< "   cenpack"
						<< std::endl;



	ScoreFunction
		sf_env,sf_pair,sf_vdw,sf_hs,sf_ss,sf_sheet,sf_cbeta,
		sf_rsigma,sf_hb_srbb,sf_hb_lrbb,sf_rg,sf_co,sf_rama,
		sf_cenpack;

	sf_env    .set_weight( scoring::env     , 1.0 );
	sf_pair   .set_weight( scoring::pair    , 1.0 );
	sf_vdw    .set_weight( scoring::vdw     , 1.0 );
	// sf_hs     .set_weight( scoring::hs      , 1.0 );
	// sf_ss     .set_weight( scoring::ss      , 1.0 );
	// sf_sheet  .set_weight( scoring::sheet   , 1.0 );
	sf_cbeta  .set_weight( scoring::cbeta   , 1.0 );
	// sf_rsigma .set_weight( scoring::rsigma  , 1.0 );
	// sf_hb_srbb.set_weight( scoring::hb_srbb , 1.0 );
	// sf_hb_lrbb.set_weight( scoring::hb_lrbb , 1.0 );
	sf_rg     .set_weight( scoring::rg      , 1.0 );
	sf_rama   .set_weight( scoring::rama    , 1.0 );
	sf_cenpack  .set_weight( scoring::cenpack   , 1.0 );

	ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	for( Size ii = 1; ii <= pdbs.size(); ++ii ) {
		string pdb( pdbs[ii] );
		Pose my_pose;
		core::io::pdb::pose_from_pdb( my_pose, *rsd_set, "centroid_pdb/"+pdb+"_0001.pdb" );

		Real s_env     = sf_env    (my_pose);
		Real s_pair    = sf_pair   (my_pose);
		Real s_vdw     = sf_vdw    (my_pose);
		Real s_hs      = sf_hs     (my_pose);
		Real s_ss      = sf_ss     (my_pose);
		Real s_sheet   = sf_sheet  (my_pose);
		Real s_cbeta   = sf_cbeta  (my_pose);
		Real s_rsigma  = sf_rsigma (my_pose);
		Real s_hb_srbb = sf_hb_srbb(my_pose);
		Real s_hb_lrbb = sf_hb_lrbb(my_pose);
		Real s_rg      = sf_rg     (my_pose);
		Real s_rama    = sf_rama   (my_pose);
		Real s_cenpack = sf_cenpack(my_pose);

		std::cout << "SCORE "
							<< pdbs[ii] << " "
							<< F( 10, 5, s_env    )
							<< F( 10, 5, s_pair   )
							<< F( 10, 5, s_vdw    )
							<< F( 10, 5, s_hs     )
							<< F( 10, 5, s_ss     )
							<< F( 10, 5, s_sheet  )
							<< F( 10, 5, s_cbeta  )
							<< F( 10, 5, s_rsigma )
							<< F( 10, 5, s_hb_srbb)
							<< F( 10, 5, s_hb_lrbb)
							<< F( 10, 5, s_rg     )
							<< F( 10, 5, s_rama   )
							<< F( 10, 5, s_cenpack   )
							<< std::endl;

	}

}


int main (int argc, char *argv[])
{

	core::init( argc, argv );

	print_scores();

	return 0;

}
