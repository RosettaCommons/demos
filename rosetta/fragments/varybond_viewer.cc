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
//#include <core/scoring/ScoringManager.hh>
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
#include "core/kinematics/Stub.hh"
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include "core/id/DOF_ID.hh"

#include <core/io/pdb/pose_io.hh>

#include <core/mm/mm_torsion_library.hh>
#include <core/mm/mm_torsion_library.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/option.hh>

#include <core/util/basic.hh>
#include <core/util/prof.hh>

#include <core/io/database/open.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Options.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include "numeric/conversions.hh"
#include "numeric/xyz.functions.hh"
#include "numeric/xyzMatrix.hh"
#include "numeric/xyzMatrix.io.hh"

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "time.h"

/////////////////////////////// TESTING
// #ifdef MAC
// #include <GLUT/glut.h>
// #else
// #include "GL/glut.h"
// #endif

// #include <pthread.h>
/////////////////////////////// TESTING
//
// typedef numeric::xyzVector<core::Real> V;
// typedef numeric::xyzMatrix<core::Real> M;

using namespace utility;
using namespace numeric;
using namespace numeric::conversions;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using core::Real;

//silly using/typedef

using namespace core;

using utility::vector1;



#include "core/fragments/Frag.hh"
#include "core/fragments/VectorFragSet.hh"
#include "core/fragments/frag_io.hh"

#include "utility/Filter.hh"
#include "utility/FilterInequality.hh"

using namespace fragments;
using namespace pose;





void
minibinitio( PoseOP p
					 , FragSetOP frags3
					 , FragSetOP frags9
					 )
{
	using namespace pose;
	using namespace protocols;
	using namespace protocols::moves;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace scoring;
	using namespace optimization;
	using core::Real;

	bool const  choose_frag_ss_check = true;
	Real const init_temp = 2.0;
	int  const increasecycles = 1;

	int const score0_cycles(  static_cast< int > (2000 * increasecycles) );
	int const score1_cycles(  static_cast< int > (2000 * increasecycles) );
	int const score25_cycles( static_cast< int > (2000 * increasecycles) );
	int const score3_cycles(  static_cast< int > (4000 * increasecycles) );

	scoring::ScoreFunction score0, score1, score2, score3, score5;

	score0.set_weight( scoring::vdw  , 1.0 );

	score1.set_weight( scoring::vdw  , 1.0 );
	score1.set_weight( scoring::env  , 1.0 );
	score1.set_weight( scoring::pair , 1.0 );

	score2.set_weight( scoring::vdw  , 1.0 );
	score2.set_weight( scoring::env  , 1.0 );
	score2.set_weight( scoring::pair , 1.0 );
	score2.set_weight( scoring::cbeta, 0.5 );

	score5.set_weight( scoring::vdw  , 1.0 );
	score5.set_weight( scoring::env  , 1.0 );
	score5.set_weight( scoring::pair , 1.0 );
	score5.set_weight( scoring::cbeta, 0.5 );

	score3.set_weight( scoring::vdw  , 1.0 );
	score3.set_weight( scoring::env  , 1.0 );
	score3.set_weight( scoring::pair , 1.0 );
	score3.set_weight( scoring::cbeta, 1.0 );

	Filter<Frag> f;

	MonteCarlo mc( p, &score0, init_temp );
	mc.set_autotemp( true, init_temp );

	std::cerr << "score0 moves" << std::endl;
	mc.score_function( &score0 );
	mc.set_temperature( init_temp );
	for( int j = 1; j <= score0_cycles; ++j ) {
		// skipping some frag picking logic like do_ss_check
		Frag f = frags9->sample();
		if( f.check(p) ) f.insert( p );
		mc.boltzmann( p );
	}

	std::cerr << "score1 moves" << std::endl;
	mc.score_function( &score1 );
	mc.set_temperature( init_temp );
	for( int j = 1; j <= score1_cycles; ++j ) {
		// skipping some frag picking logic like do_ss_check
		// CHECKPOINT!!!!
		Frag f = frags9->sample();
		if( f.check(p) ) f.insert( p );
		mc.boltzmann( p );
	}

	std::cerr << "score25 moves" << std::endl;
	mc.set_temperature( init_temp );
	for( int jj = 1; jj <= 1; ++jj ) {
		for( int kk = 1; kk <= 10; ++kk ) {
			if( kk%2 == 0 || kk > 7 ) { // score 2
				// chainbreak = jj*kk*0.25
				mc.score_function(&score2);
			} else { // score 5
				// chainbreak = jj*kk*0.05
				mc.score_function(&score5);
			}
			for( int j = 1; j <= score25_cycles; ++j ) {
				// skipping some frag picking logic like do_ss_check
				// CHECKPOINT!!!!
				Frag f = frags9->sample();
				if( f.check(p) ) f.insert( p );
			}
		}
	}

	std::cerr << "score3 moves" << std::endl;
	mc.score_function( &score3 );
	for( int kk = 1; kk <= 3; ++kk ) {
		if( 2 == kk ) ; // choose_frag_set_top_N_frags(number3merfrags);
		for( int j = 1; j <= score3_cycles; ++j ) {
			// CHECKPOINT
			if( 1 == kk ) {
				Frag f = frags3->sample();
				if( f.check(p) ) f.insert( p );
			} else {
				; // choose_fragment_gunn_pose( pose, size, gunn_cutoff );
			}
			mc.boltzmann( p );
		}
	}

}

struct FragSize {
	int get_numeric( Frag f ) { return f.fragData()->size(); }
	std::string description() { return "max element of vector"; }
};


void*
my_main( void* ) {

	//assert( test_Transform(99999) );

	// frag_test();
	//   exit(-1);

	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace scoring;
	using namespace optimization;
	using namespace core::io::pdb;
	using utility::Filter;

	ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	pose::Pose cenpose;
	core::io::pdb::pose_from_pdb( cenpose, *rsd_set, "input/1ubi.pdb" );

	protocols::viewer::add_conformation_viewer( cenpose.conformation(), "start_pose" );

	FragSetOP frags( new VectorFragSet ), frags_ideal( new VectorFragSet );
	std::cout << "reading frags" << std::endl;
	add_frags_from_file(frags,"input/1ubi.will_frags");
	std::cout << "reading ideal frags" << std::endl;
	add_frags_from_file(frags_ideal,"input/1ubi_ideal.will_frags");

	FilterInequalityGenerator<Frag,int,FragSize> len;
	FragSetOP frags3 = frags_ideal->subset( len == 3 );
	FragSetOP frags9 = frags_ideal->subset( len == 9 );

	for( int ii = 1; ii <= cenpose.total_residue(); ++ii ) {
		cenpose.set_phi(ii, 180 );
		cenpose.set_psi(ii, 180 );
	}
	minibinitio( &cenpose, frags3, frags9 );

	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	// void * tmp;
	// my_main(tmp);
}









