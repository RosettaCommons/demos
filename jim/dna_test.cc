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
// some silly helper routines:
//




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
//
// if position i in seq1 is aligned with position j in seq2, mapping[ i ] == j
// if position i in seq1 is unaligned, mapping[ i ] == 0
//
void
read_alignment_file(
	std::string const & filename,
	std::string & seq1,
	std::string & seq2,
	utility::vector1< int > & mapping // from numbering in sequence 1 to numbering in sequence 2
)
{
	std::string align1, align2;
	{ // parse the file
		std::ifstream data( filename.c_str() );
		std::string line;
		// 1st sequence
		getline( data,line );
		assert( line[0] == '>' );
		getline( data, align1 );
		// 2nd sequence
		getline( data, line );
		assert( line[0] == '>' );
		getline( data, align2 );
		data.close();
	}

	assert( align1.size() == align2.size() );

	seq1.clear();
	seq2.clear();
	mapping.clear();
	int pos1(0), pos2(0);
	for ( Size i=1; i<= align1.size(); ++i ) {
		char const al1( align1[i] ), al2( align2[i] );
		bool const gap1( al1 == '.' || al1 == '-' );
		bool const gap2( al2 == '.' || al2 == '-' );
		if ( !gap2 ) {
			++pos2;
			seq2 += al2;
		}
		if ( !gap1 ) {
			++pos1;
			seq1 += al1;
			if ( !gap2 ) {
				mapping.push_back( pos2 );
			} else {
				mapping.push_back( 0 ); // unaligned
			}
		}
	}

	assert( mapping.size() == seq1.size() );
}

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

	Real const CUTOFF( 4.5 );

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

		Pose pose;

		pose_from_pdb( pose, filename );

		std::cout << "read file: " << filename << ' '<< pose.total_residue() << std::endl;

		std::cout << "SEQUENCE: " << filename << ' ' << pose.sequence() << std::endl;

		for ( Size i=1; i<= pose.total_residue(); ++i ) {

			Residue const & rsd1( pose.residue(i) );
			if ( rsd1.is_protein() ) {

				bool is_interface( false );

				for ( Size j=1; j<= pose.total_residue() && !is_interface; ++j ) {

					Residue const & rsd2( pose.residue(j) );
					if ( rsd2.is_DNA() ) {

						for ( Size ii=1; ii<= rsd1.natoms() && !is_interface; ++ii ) {
							for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
								Vector const & xyz1( rsd1.xyz( ii ) );
								Vector const & xyz2( rsd2.xyz( jj ) );
								if ( xyz1.distance( xyz2 ) < CUTOFF ) {
									is_interface = true;
									break;
								}
							}
						}
					}
				}

				if ( is_interface ) std::cout << "interface residue: " << filename << ' ' << i << ' ' <<
															rsd1.aa() << std::endl;
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
void
simple_dna_test()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;

	Pose pose;


	ResidueTypeSet const & residue_set( *(ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) ) );
	pose_from_pdb( pose, residue_set, "input/1aay.pdb" );

	std::cout << pose.sequence() << std::endl;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 ); // LJ attractive
	scorefxn.set_weight( fa_rep, 0.44 ); // LJ repulsize
	scorefxn.set_weight( fa_sol, 0.65 ); // LK solvation
	//scorefxn.set_weight( fa_pair, 0.49 ); // knowledge-based rsd-rsd electrostatics
	scorefxn.set_weight( fa_dun, 0.56 ); // Dunbrack rotamer energy
	scorefxn.set_weight( rama, 0.2 ); // ramachandran score
	scorefxn.set_weight( hbond_lr_bb, 1.17 ); // long-range backbone hbonds
	scorefxn.set_weight( hbond_sr_bb, 1.17 ); // short-range (helical) backbone hbonds
	scorefxn.set_weight( hbond_bb_sc, 1.17 ); // backbone-sidechain hbonds
	scorefxn.set_weight( hbond_sc   , 1.10 ); // sidechain-sidechain hbonds

	std::cout << "score1: " << scorefxn( pose ) << std::endl;


	{ // try rotamer optimizations //////////////////////////////////////////////////////////////////////////////
		Energy score_orig = scorefxn( pose );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		// new interface for packer task: task starts out redesigning at all protein positions with all amino
		// acids, and starts out repacking at all other residues -- the new interface requires the protocol
		// make restrictions from here, providing commutativity of restrictions.
		task->initialize_from_command_line();
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		pack::rotamer_trials( pose, scorefxn, task );

		Energy rottrial_score = scorefxn( pose );

		std::cout << "Completed rotamer_trials_test() with new score: " << rottrial_score << " vs orig: " <<
			score_orig << std::endl;

		dump_pdb( pose, "test_rottrials.pdb" );

		// now try packing
		task->or_include_current( true ); // tmp hack
		pack::pack_rotamers( pose, scorefxn, task);
		Energy pack_score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() with new score: " << pack_score << " vs orig: " <<
			rottrial_score << std::endl;

		dump_pdb( pose, "test_packrots.pdb" );

		// these are useful but cost a little time to get
		scorefxn.accumulate_residue_total_energies( pose );

		pose.energies().show( std::cout );

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
