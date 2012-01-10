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

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>


#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/util.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

using namespace core;

using utility::vector1;


typedef std::map< std::string, Vector > ResidueCoords;
typedef std::map< std::string, ResidueCoords > Coords;

typedef vector1< std::string > Strings;

///////////////////////////////////////////////////////////////////////////////
std::ostream &
operator<< ( std::ostream & out, Vector const & v ) {
	out << "( " << F(9,3,v(1)) << " , " << F(9,3,v(2)) << " , " << F(9,3,v(3)) << " )";
	return out;
}




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void
simple_dna_test2()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;


	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace optimization;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );

	scorefxn.set_weight( h2o_intra, 0.01 );
	scorefxn.set_weight( h2o_hbond, 1.0 );

	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

//	scorefxn.set_weight( dna_bs, 0.5 );
//	scorefxn.set_weight( dna_bp, 0.5 );

	// read list file of pdbs

	utility::vector1< std::string > filenames;

	{
		std::string const listfile( start_file() );
		std::ifstream data( listfile.c_str() );
		std::string line;
		while ( getline( data,line ) ){
			filenames.push_back( line );
		}
		data.close();
	}

	for ( Size n=1; n<= filenames.size(); ++n ) {
		std::string const filename( filenames[n] );

		std::cout << "Working on file: " << filename << std::endl;

		Pose pose;

		pose_from_pdb( pose, filename );

		std::cout << "read file: " << filename << ' '<< pose.total_residue() << std::endl;

		std::cout << "SEQUENCE: " << filename << ' ' << pose.sequence() << std::endl;

		Energy score_orig = scorefxn( pose );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
		task->set_bump_check( true );

		clock_t starttime = clock();
		pack::pack_rotamers( pose, scorefxn, task);
		clock_t stoptime = clock();
		std::cout << "pack_rotamers took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
		Energy score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() #1 with new score: " << score << " vs orig: " << score_orig << std::endl;

		dump_pdb( pose, "post_packrots.pdb" );


		// Test minimization with adduct waters
		//
		kinematics::MoveMap mm;
		mm.set_bb( true );
		mm.set_chi( true );

		{
			// setup the options
			MinimizerOptions options( "dfpmin", 0.01, true /*use_nblist*/, true /*deriv_check*/ );

			AtomTreeMinimizer minimizer;
			minimizer.run( pose, mm, scorefxn, options );
			dump_pdb( pose, "post_minimization.pdb" );
		}

		scorefxn.show( std::cout );

	}
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	//using namespace core;
	core::init( argc, argv );

	simple_dna_test2();

}
