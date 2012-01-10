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

#define SHEFFLERDEBUG


//#include "boost/shared_ptr.hpp"


// libRosetta headers
//#include <core/options/option.hh>

#include <core/types.hh>

#include <core/chemical/AtomSet.hh>
#include <core/chemical/MMAtomSet.hh>

#include <core/conformation/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueSet.hh>
#include <core/conformation/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/residue_io.hh>

#include <core/scoring/Etable.hh>
#include <core/scoring/KnowledgeManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/count_pair/CountPairFunction.hh>
#include <core/scoring/methods/EnergyMethodManager.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomID_Map.hh>
#include <core/kinematics/AtomID_Map.Pose.hh>

#include <core/mm/mm_torsion_library.hh>
#include <core/mm/mm_torsion_library.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/option.hh>

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


#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;



using namespace core;

using utility::vector1;


typedef std::map< std::string, std::map< std::string, numeric::xyzVector< Real > > > Coords;

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

void
dump_pdb(
	pose::Pose const & pose,
	std::ostream & out
)
{
	//using namespace core;
	int const nres( pose.total_residue() );

	int number(0);


	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );
			++number;
			out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << "  " << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(6,2,1.0) << F(6,2,1.0) << '\n';
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
dump_pdb(
	pose::Pose const & pose,
	std::string const & filename
)
{
	std::ofstream out( filename.c_str() );
	dump_pdb( pose, out );
	out.close();
}



///////////////////////////////////////////////////////////////////////////////
// super-simple pdb reader
//

void
read_pdb(
	std::string const & filename,
	Strings & resids,
	Strings & sequence,
	Coords & coords
	)
{
	using ObjexxFCL::float_of;
	typedef numeric::xyzVector< Real > Vector;

	resids.clear();
	sequence.clear();
	coords.clear();


	std::ifstream data( filename.c_str() );
	std::string line;
	while ( getline( data, line ) ) {
		if ( line.substr(0,6) == "ATOM  " ) {
			// parse the info
			std::string const atom_name( line.substr(12,4));
			std::string const name3( line.substr(17,3));
			std::string const resid( line.substr(22,5));

			Vector const xyz
				( float_of( line.substr(30,8) ),
					float_of( line.substr(38,8) ),
					float_of( line.substr(46,8) ) );

			if ( coords.find( resid ) == coords.end() ) {
				resids.push_back( resid );
				sequence.push_back( name3 );
			}

			coords[ resid ] [ atom_name ] = xyz;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
// super-simple
//
void
pose_from_pdb(
							pose::Pose & pose,
							conformation::ResidueSet const & residue_set,
							std::string const & filename
							)
{
	//using namespace core;
	using namespace conformation;

	typedef numeric::xyzVector< Real > Vector;

	// reset current data
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:222 (" << ")" << std::endl;
//#endif
	pose.clear();

	Coords coords;
	Strings resids, sequence;

//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:227 (" << ")" << std::endl;
//#endif
	read_pdb( filename, resids, sequence, coords );
	int const nres_pdb( resids.size() );
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:232 (" << ")" << std::endl;
//#endif
	for ( int i=1; i<= nres_pdb; ++i ) {
		std::string const pdb_name( sequence[i] );
		std::string const resid( resids[i] );
		std::map< std::string, Vector > const & xyz
			( coords.find( resid )->second );

		// which residues match this 3letter name?
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:242 (" << ")" << std::endl;
//#endif
		ResidueTypeCAPs const & rsd_type_list( residue_set.name3_map( pdb_name ) );
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:249 (" << ")" << std::endl;
//#endif
		if ( rsd_type_list.empty() ) {
			std::cout << "Unrecognized aa: " << pdb_name << '\n';
			continue;
		}
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:256 (" << ")" << std::endl;
//#endif

		// look for perfect match:
		bool matched( false );
		for ( Size j=1; j<= rsd_type_list.size(); ++j ) {
			ResidueType const & rsd_type( *(rsd_type_list[j]) );

			if ( Size( rsd_type.natoms() ) != xyz.size() ) continue;
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:260 (" << ")" << std::endl;
//#endif
			int mismatches(0);
			for ( Size k=1; k<= xyz.size(); ++k ) {
				if ( xyz.count( rsd_type.atom_name(k) ) == 0 ) ++mismatches;
			}
			if ( mismatches ) continue;

			matched = true;

			// found a perfect match! fill in the coords
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:272 (" << ")" << std::endl;
//#endif
			ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );

			for ( Size k=1; k<= xyz.size(); ++k ) {
				new_rsd->atom(k).xyz( xyz.find( new_rsd->atom_name(k) )->second );
			}

			pose.append_residue( new_rsd );
		} // j=1,rsd_type_list.size()
		if ( !matched ) {
			std::cout << "Unrecognized residue: " << pdb_name << std::endl;
		}

	} // i=1,nres_pdb
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:288 (" << ")" << std::endl;
//#endif

	pose.fold_tree( kinematics::FoldTree( pose.total_residue() ) );
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:293 (" << ")" << std::endl;
//#endif

}


///////////////////////////////////////////////////////////////////////////////
void
setup_atom_number(
	pose::Pose const & pose,
	kinematics::AtomID_Map< int > & atom_number
)
{
	kinematics::initialize( atom_number, pose );

	int number(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		for ( Size j=1; j<= Size(pose.residue(i).natoms()); ++j ) {
			atom_number[ kinematics::AtomID( j, i ) ] = ++number;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
show_rasmol_hbonds(
	scoring::hbonds::HBondSet const & hbond_set,
	kinematics::AtomID_Map< int > const & atom_number,
	std::ostream & out
)
{
	using namespace scoring::hbonds;

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		if ( hbond_set.allow_hbond(i) ) {
			HBond const & hb( hbond_set.hbond(i) );
			kinematics::AtomID const hatm( hb.don_hatm(), hb.don_res() );
			kinematics::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
			out << "monitor " << atom_number[ hatm ] << ' ' << atom_number[ aatm ] <<
				'\n';
		}
	}

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
test_dunbrack_io(
	conformation::ResidueSet const & residue_set
)
{
	//using namespace core;
	using namespace scoring;
	using namespace scoring::dunbrack;

	RotamerLibrary const & rot_lib
		( KnowledgeManager::get_instance()->get_RotamerLibrary() );

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		rot_lib.rotamer_energy( pose.residue(i) );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
test_rama(
	conformation::ResidueSet const & residue_set
)
{
	using namespace scoring;
	std::cout << "test rama" << std::endl;

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	ScoreFunction scorefxn;
	scorefxn.set_weight( rama , 0.2 );

	scorefxn( pose );

	std::cout << "rama succeeded" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
simple_rotamer_test(
	conformation::ResidueSet const & residue_set
)
{
	using namespace conformation;

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	ResidueOP ile_rsd( pose.residue(3).clone() );
	ResidueOP pro_rsd( pose.residue(70).clone() );

	assert( pro_rsd->name() == "PRO" && ile_rsd->name() == "ILE" );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.nchi() >= 1 ) {
			Residue rot( rsd.type(), rsd, pose.conformation() );
			rot.set_chi( 1, 60.0 );

			Residue rot2( rsd.type(), *ile_rsd, pose.conformation() );
			rot2.set_chi( 1, 60.0 );

			// this fails b/c PRO has no backbone H
			Residue rot3( rsd.type(), *pro_rsd, pose.conformation() );
			rot3.set_chi( 1, 60.0 );
		}

		if ( true ) { //rsd.name() != "PRO" && rsd.name() != "GLY" ) {
			Residue rot( ile_rsd->type(), rsd, pose.conformation() );
			rot.set_chi(1,180.0);
			rot.set_chi(2,180.0);
			ResidueOP new_rsd( rot.create_residue() );

			pose.replace_residue( i, *new_rsd, false );

		}
	}

	dump_pdb( pose, "test_out.pdb" );

	//exit(0);
}
///////////////////////////////////////////////////////////////////////////////
void
simple_conformation_test(
	conformation::ResidueSet const & residue_set
)
{


	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	pose.set_phi( 3, pose.phi(3) );

	dump_pdb( pose, "output/test1.pdb" );

	pose.set_phi( 3, 0.0 );

	dump_pdb( pose, "output/test2.pdb" );

	pose.set_chi( 1, 3, 60.0 );

	dump_pdb( pose, "output/test3.pdb" );

	Size const nres( pose.total_residue() );
	kinematics::FoldTree f( nres );
	f.reorder( nres );

	pose.fold_tree( f );

	pose.set_phi( 3, 120.0 );

	dump_pdb( pose, "output/test4.pdb" );


}

///////////////////////////////////////////////////////////////////////////////

void
simple_dna_test(
	conformation::ResidueSet const & residue_set
)
{

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in_dna.pdb" );

	// setup scorefxn
	scoring::ScoreFunction scorefxn;
	scorefxn.set_weight( scoring::fa_atr, 0.80 );
	scorefxn.set_weight( scoring::fa_rep, 0.65 );
	scorefxn.set_weight( scoring::fa_sol, 0.44 );

	scorefxn( pose );

	dump_pdb(pose, "junk_dna.pdb" );
}

///////////////////////////////////////////////////////////////////////////////
void
simple_min_test(
	conformation::ResidueSet const & residue_set
)
{
	using namespace optimization;
	using pose::Pose;
	using kinematics::AtomID;
	using kinematics::DOF_ID;
	using kinematics::PHI;
	using kinematics::THETA;
	using kinematics::D;

	pose::Pose start_pose;
	pose_from_pdb( start_pose, residue_set, "test_in.pdb" );


	{ // test out a couple different scoring functions
		kinematics::MoveMap mm;

		// setup moving dofs
		for ( int i=30; i<= 35; ++i ) {
			mm.set_bb ( i, true );
			mm.set_chi( i, true );
		}

		// setup the options
		MinimizerOptions options( "linmin", 10.0, true /*use_nblist*/,
			true /*deriv_check*/ );

		AtomTreeMinimizer minimizer;
		{ // just fa_rama
			scoring::ScoreFunction scorefxn;
			scorefxn.set_weight( scoring::rama, 1.0 );

			Pose pose;
			pose = start_pose;
			std::cout << "MINTEST: rama" << std::endl;
			minimizer.run( pose, mm, scorefxn, options );
		}

		/////////////
		return;
		/////////////

		{ // just fa_dun
			scoring::ScoreFunction scorefxn;
			scorefxn.set_weight( scoring::fa_dun, 1.0 );

			Pose pose;
			pose = start_pose;
			std::cout << "MINTEST: fa_dun" << std::endl;
			minimizer.run( pose, mm, scorefxn, options );
		}

		{ // just fa_atr
			scoring::ScoreFunction scorefxn;
			scorefxn.set_weight( scoring::fa_atr, 0.80 );

			Pose pose;
			pose = start_pose;
			std::cout << "MINTEST: atr" << std::endl;
			minimizer.run( pose, mm, scorefxn, options );
		}

		{ // just fa_atr, rep, sol
			scoring::ScoreFunction scorefxn;
			scorefxn.set_weight( scoring::fa_atr, 0.80 );
			scorefxn.set_weight( scoring::fa_rep, 0.44 );
			scorefxn.set_weight( scoring::fa_sol, 0.65 );

			Pose pose;
			pose = start_pose;
			std::cout << "MINTEST: atr-rep-sol" << std::endl;
			minimizer.run( pose, mm, scorefxn, options );
		}

		{ // just backbone hbonds
			scoring::ScoreFunction scorefxn;
			scorefxn.set_weight( scoring::hbond_lr_bb, 1.0 );
			scorefxn.set_weight( scoring::hbond_sr_bb, 1.0 );

			Pose pose;
			pose = start_pose;
			std::cout << "MINTEST: bb hbonds" << std::endl;
			minimizer.run( pose, mm, scorefxn, options );
		}

		{ // all hbonds
			scoring::ScoreFunction scorefxn;
			scorefxn.set_weight( scoring::hbond_lr_bb, 1.0 );
			scorefxn.set_weight( scoring::hbond_sr_bb, 1.0 );
			scorefxn.set_weight( scoring::hbond_bb_sc, 1.0 );
			scorefxn.set_weight( scoring::hbond_sc, 1.0 );

			Pose pose;
			pose = start_pose;
			std::cout << "MINTEST: all hbonds" << std::endl;
			minimizer.run( pose, mm, scorefxn, options );
		}

		return;

	} // scope


	Pose pose;
	pose = start_pose;

	// set the moving dofs
	kinematics::MoveMap mm1, mm2, mm3, mm4;
	// single backbone
	mm1.set_bb( 4, true );

	// all bb and chi
	mm2.set_bb( true );
	mm2.set_chi( true );

	// single dof
	mm3.set( DOF_ID( AtomID(1,4), PHI ), true );

	// everything!
	mm4.set( PHI, true );
	mm4.set( THETA, true );
	mm4.set( D, true );

	// setup scorefxn
	scoring::ScoreFunction scorefxn;
	scorefxn.set_weight( scoring::fa_atr, 0.80 );

	MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/,
														true /*deriv_check*/ );

	AtomTreeMinimizer minimizer;

	dump_pdb( pose, "output/before.pdb" );

	minimizer.run( pose, mm1, scorefxn, options );
	dump_pdb( pose, "output/after1.pdb" );

	exit(0);

	minimizer.run( pose, mm2, scorefxn, options );
	dump_pdb( pose, "output/after2.pdb" );

	minimizer.run( pose, mm3, scorefxn, options );
	dump_pdb( pose, "output/after3.pdb" );

	minimizer.run( pose, mm4, scorefxn, options );
	dump_pdb( pose, "output/after4.pdb" );

}
///////////////////////////////////////////////////////////////////////////////
void
dump_hbond_pdb(
	pose::Pose & pose,
	std::string const & filename
)
{
	using namespace optimization;
	using kinematics::AtomID;
	using kinematics::DOF_ID;
	using kinematics::PHI;
	using kinematics::THETA;
	using kinematics::D;

	//
	scoring::hbonds::HBondSet hbond_set;
	pose.update_residue_neighbors();
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	// setup the global atom numbering that would be used for pdb output
	kinematics::AtomID_Map< int > atom_number;
	setup_atom_number( pose, atom_number );

	std::ofstream out( filename.c_str() );
	out << "load pdb inline\n";
	show_rasmol_hbonds( hbond_set, atom_number, out );
	out << "exit\n";
	dump_pdb( pose, out );
	out.close();
}


///////////////////////////////////////////////////////////////////////////////
void
simple_copy_test(
	conformation::ResidueSet const & residue_set
)
{

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	pose::Pose new_pose;
	new_pose = pose;

	dump_pdb( new_pose, "test_copy.pdb" );
}


///////////////////////////////////////////////////////////////////////////////
void
simple_hbond_test(
	conformation::ResidueSet const & residue_set
)
{
	using namespace optimization;
	using kinematics::AtomID;
	using kinematics::DOF_ID;
	using kinematics::PHI;
	using kinematics::THETA;
	using kinematics::D;

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	dump_hbond_pdb( pose, "test_out.hbonds.pdb" );

	// setup scorefxn
	scoring::ScoreFunction scorefxn;
	scorefxn.set_weight( scoring::hbond_lr_bb, 1.0 );
	scorefxn.set_weight( scoring::hbond_sr_bb, 1.0 );
	scorefxn.set_weight( scoring::hbond_bb_sc, 1.0 );
	scorefxn.set_weight( scoring::hbond_sc, 1.0 );

	scorefxn( pose );

	//pose.energies().show( std::cout );

	// set the moving dofs
	kinematics::MoveMap mm1, mm2, mm3, mm4;

	// single backbone
	mm1.set_bb( 4, true );

	// all bb and chi
	mm2.set_bb( true );
	mm2.set_chi( true );

	// single dof
	mm3.set( DOF_ID( AtomID(1,4), PHI ), true );

	// everything!
	mm4.set( PHI, true );
	mm4.set( THETA, true );
	mm4.set( D, true );

	MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/,
														true /*deriv_check*/ );

	AtomTreeMinimizer minimizer;

	dump_pdb( pose, "output/before.pdb" );

	//minimizer.run( pose, mm1, scorefxn, options );
	//dump_pdb( pose, "output/after1.pdb" );

	minimizer.run( pose, mm2, scorefxn, options );
	dump_pdb( pose, "output/after2.pdb" );

	minimizer.run( pose, mm3, scorefxn, options );
	dump_pdb( pose, "output/after3.pdb" );

	minimizer.run( pose, mm4, scorefxn, options );
	dump_pdb( pose, "output/after4.pdb" );

}


///////////////////////////////////////////////////////////////////////////////
void
atom_tree_torsion_test(
	conformation::ResidueSet const & residue_set
)
{
	using kinematics::AtomID;

	pose::Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	pose.conformation().debug_residue_torsions();

	kinematics::AtomTree const & atom_tree( pose.atom_tree() );

	for ( Size i=2; i<= pose.total_residue() - 1; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		assert( Size( rsd.seqpos() ) == i );
		for ( Size chino=1; chino<= rsd.nchi(); ++chino ) {
			AtomID atom2( rsd.chi_atoms()[chino][2], i );
			AtomID atom3( rsd.chi_atoms()[chino][3], i );

			// loop over atoms bonded to atom2
			vector1< int > atom2_nbrs( rsd.bonded_neighbor( atom2.atomno() ) );
			vector1< int > atom3_nbrs( rsd.bonded_neighbor( atom3.atomno() ) );
			for ( Size ii=1; ii<= atom2_nbrs.size(); ++ii ) {
				for ( Size jj=1; jj<= atom3_nbrs.size(); ++jj ) {
					AtomID const atom1( atom2_nbrs[ii], i );
					AtomID const atom4( atom3_nbrs[jj], i );
					if ( atom1 == atom3 || atom2 == atom4 ) continue;
					Real offset;
					atom_tree.torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset);
				}
			}
		}
	}


	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );

		for ( int r=1; r<= 2; ++r ) {
			kinematics::TorsionType const type( r == 1 ? kinematics::BB :
																					kinematics::CHI );
			Size const n( r == 1 ? rsd.mainchain_atoms().size() : rsd.nchi() );

			for ( Size j=1; j<= n; ++j ) {
				kinematics::TorsionID const tor_id( i, type, j );
				//std::cout << "set_torsion: " << tor_id <<
				//	" =============================================" << std::endl;
				pose.set_torsion( tor_id, 180.0 );
				pose.conformation().debug_residue_torsions();
				assert( std::abs( util::subtract_degree_angles( pose.torsion( tor_id ),
																												180.0 ) ) < 1e-3 );
			}
		}
	}
	dump_pdb( pose, "extended.pdb" );
}

///////////////////////////////////////////////////////////////////////////////
void
fa_scorefxn_test(
	conformation::ResidueSet const & residue_set
)
{
	// setup etable, fa scorefxn
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;



	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );
	scorefxn.set_weight( fa_dun, 0.49 );
	scorefxn.set_weight( rama, 0.2 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

	Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );

	std::cout << "pose-score: " << scorefxn( pose ) << std::endl;

	scorefxn.accumulate_residue_total_energies( pose );

}

///////////////////////////////////////////////////////////////////////////////
void
rotamer_trials_test(
	conformation::ResidueSet const & residue_set
)
{
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );

	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

	std::cout << "atr-weight: " << scorefxn[fa_atr] << std::endl;
	std::cout << "rep-weight: " << scorefxn[fa_rep] << std::endl;


	Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );
	Energy score_orig = scorefxn( pose );

	//for ( int jj = 1; jj <= pose.total_residue(); ++jj )
	//{
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose.total_residue() ));
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii )
	{
		task->set_pack_residue( ii, true );
	}
	//task->set_pack_residue( jj , true );

	pack::rotamer_trials( pose, scorefxn, task, residue_set );

	Energy score = scorefxn( pose );

	std::cout << "Completed rotamer_trials_test() with new score: " << score << " vs orig: " << score_orig << std::endl;

	pose.energies().show( std::cout );

	//std::string outname = "test_rotrials_";
	//std::stringstream ss;
	//ss << jj;
	//outname += ss.str();
	//outname += ".pdb";

	//std::cout << outname << std::endl;

	//dump_pdb( pose, outname.c_str() );
	dump_pdb( pose, "rotamer_trials_all.pdb");
	//}

}

///////////////////////////////////////////////////////////////////////////////
void
pack_rotamers_test(
	conformation::ResidueSet const & residue_set
)
{
	//using namespace core;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	ScoreFunction scorefxn;
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );

	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );

	Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );
	Energy score_orig = scorefxn( pose );

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose.total_residue() ));
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii )
	{
		task->set_pack_residue( ii, true );
	}
	task->set_include_current( true );

	pack::pack_rotamers( pose, scorefxn, residue_set, task);
	Energy score = scorefxn( pose );

	std::cout << "Completed pack_rotamers_test() with new score: " << score << " vs orig: " << score_orig << std::endl;

	dump_pdb( pose, "test_packrots.pdb" );

}

///////////////////////////////////////////////////////////////////////////////
void
patch_test(
	conformation::ResidueSet & residue_set
) {
	using namespace conformation;
	using namespace pose;

	residue_set.apply_patches( io::database::full_name( "patches.txt" ) );

	Pose pose;
	pose_from_pdb( pose, residue_set, "test_in.pdb" );
	Size const nres( pose.total_residue() );

	if ( false ) {
	for ( Size i=2; i<= nres-1; ++i ) {
		Residue const & rsd( pose.residue(i) );

		ResidueOP new_rsd( rsd.clone() );

		for ( Size j=1; j<= rsd.natoms(); ++j ) {

			bool is_chi_atom( false );
			for ( Size k=1; k<= rsd.nchi(); ++k ) {
				if ( j == rsd.chi_atoms(k)[4] ) is_chi_atom = true;
			}
			if ( is_chi_atom ) continue;

			Vector v( rsd.build_atom_ideal( j, pose.conformation() ) ), old_v( rsd.atom(j).xyz() );
			Real const dev( v.distance( old_v ) );
			if ( dev > 0.5 ) {
				std::cout << i << ' ' << rsd.name() << ' ' << j << ' ' << rsd.atom_name(j) << ' ' <<
					F(9,3,v.distance(old_v)) << ' ' << old_v << ' ' << v << std::endl;
				new_rsd->atom(j).xyz( v );
			}
		}

		//pose.replace_residue( i, *new_rsd, false );
	}
	}

	// build termini
	ResidueSelector sel1, sel2;
	sel1.add_line( "AA ASP" );
	sel1.add_line( "VARIANT_TYPE LOWER_TERMINUS" );
	sel1.add_line( "NOT VARIANT_TYPE UPPER_TERMINUS" );

	sel2.add_line( "AA LEU" );
	sel2.add_line( "PROPERTY UPPER_TERMINUS"); // alternative to variant_type in this case
	sel2.add_line( "NOT PROPERTY LOWER_TERMINUS");

	ResidueTypeCAPs m1,m2;
	residue_set.select_residues( sel1, m1 );
	residue_set.select_residues( sel2, m2 );
	std::cout << "sel1: " << m1.size() << ' ' << m1[1]->name() << std::endl;
	std::cout << "sel2: " << m2.size() << ' ' << m2[1]->name() << std::endl;

	ResidueOP nterm_asp( ResidueFactory::create_residue( *(m1[1]), pose.residue(1), pose.conformation() ) );
	ResidueOP cterm_leu( ResidueFactory::create_residue( *(m2[1]), pose.residue(nres), pose.conformation() ) );

	pose.replace_residue( 1, *nterm_asp, false );
	pose.replace_residue( nres, *cterm_leu, false );
	dump_pdb( pose, "tmp.pdb" );

	kinematics::AtomID_Mask missing;
	kinematics::initialize( missing, pose );
	for ( Size i=1; i<= nres; ++i ) {
		Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			if ( rsd.atom_type(j).is_hydrogen() ) {
				missing[ kinematics::AtomID(j,i) ] = true;
			}
		}
	}

	pose.conformation().fill_missing_atoms( missing );
	dump_pdb( pose, "tmp2.pdb" );
}


///////////////////////////////////////////////////////////////////////////////
// should run to completion
//
void simple_benchmark(	conformation::ResidueSet & residue_set ){

	simple_min_test( residue_set );

	// simple_rotamer_test( residue_set );
	//
	// rotamer_trials_test( residue_set );
	//
	// pack_rotamers_test( residue_set );
	//
	// fa_scorefxn_test( residue_set );
	//
	// test_rama( residue_set );
	//
	// simple_conformation_test( residue_set );
	//
	// atom_tree_torsion_test( residue_set );
	//
	// // this will modify the residue set
	// patch_test( residue_set );

}

void simple_benchmark_smart_ptr( conformation::ResidueSetOP rs ){
	simple_benchmark( *rs );
}
void simple_benchmark_ptr( conformation::ResidueSet *rs ){
	simple_benchmark( *rs );
}
void simple_benchmark_ref( conformation::ResidueSet & rs ){
	simple_benchmark( rs );
}
void simple_benchmark_by_value( conformation::ResidueSet rs ){
	simple_benchmark( rs );
}

void cpp_init()
{
	int argc = 3;
	char ** argv = new char*[3];
	argv[0] = new char[999];
	argv[1] = new char[999];
	argv[2] = new char[999];

	std::strcpy(argv[0],"dummy");
	std::strcpy(argv[1],"-database");
	std::strcpy(argv[2],"/Users/sheffler/svn/branches/minirosetta_database");

	std::cout << "MAKING ARGV:" << std::endl;
	std::cout << "argv[0]: '" << argv[0] << "'" << std::endl;
	std::cout << "argv[1]: '" << argv[1] << "'" << std::endl;
	std::cout << "argv[2]: '" << argv[2] << "'" << std::endl;

	options::initialize().load( argc, argv, true );
	options::process();

	numeric::random::RandomGenerator::initializeRandomGenerators(
		 111111, numeric::random::_RND_TestRun_, "ran3");

}

conformation::ResidueSetOP
res_set_ref_to_op(conformation::ResidueSet & rs) {
	return conformation::ResidueSetOP( &rs );
}

conformation::ResidueSet &
res_set_op_to_ref( conformation::ResidueSetOP rs) {
	return *rs;
}

conformation::ResidueSetOP
make_res_set() {
	chemical::AtomSetOP atom_set = new chemical::AtomSet();
	atom_set->read_file( io::database::full_name( "atom_properties.txt" ) );

	// read some MM atom properties
	chemical::MMAtomSetOP mm_atom_set = new chemical::MMAtomSet();
	mm_atom_set->read_file( io::database::full_name( "mm_atom_properties.txt" ) );

	// read in a list of MM torsion params
	mm::MMTorsionLibrary mmtl( io::database::full_name( "mm_torsion_params.txt" ), mm_atom_set);
	// std::cout << "Print all" << std::endl;
	// mmtl.pretty_print();
	// std::cout << "Print all for key: 10 10 10 10" << std::endl;
	// mmtl.pretty_print(10,10,10,10);
	// std::cout << "Print all for key CT1 CT1 CT1 CT1" << std::endl;
	// mmtl.pretty_print("CT1","CT1","CT1","CT1");

	// read a list of residues
	conformation::ResidueSetOP residue_set( new conformation::ResidueSet( atom_set, mm_atom_set ) );
	residue_set->read_list_of_residues( io::database::full_name( "rsd_list.txt" ) );

	// setup an etable
	scoring::EtableOP etable_ptr
		( new scoring::Etable( atom_set, scoring::EtableOptions() ) );

	// store it in the "global" single instance of the energy fxn manager
	scoring::methods::EnergyMethodManager::get_instance()->add_etable
		( "default", etable_ptr );

	return residue_set;
}

void
cpp_run_tests() {
	cpp_init();
	conformation::ResidueSet & rs( *make_res_set() );
	simple_benchmark( rs );
	std::cerr << "DONE!" << std::endl;
}

#define PYTHON
#ifdef PYTHON

#include "boost/python.hpp"
#include "boost/python/suite/indexing/indexing_suite.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "boost/python/suite/indexing/map_indexing_suite.hpp"

BOOST_PYTHON_MODULE_INIT(_test) {
	using namespace boost::python;
	using namespace core;
	using namespace chemical;
	using namespace conformation;
	using namespace scoring;
	using namespace methods;
	using namespace optimization;
	using namespace pose;
	using namespace kinematics;
	using std::string;

	def("options_initialize",options::initialize,return_value_policy<reference_existing_object>());
	def("options_process"   ,options::process,return_value_policy<reference_existing_object>());

	class_< std::vector<std::string> > ("VecStr")
			.def (vector_indexing_suite< std::vector<std::string> >());
	class_< std::vector<int> > ("VecInt")
			.def (vector_indexing_suite< std::vector<int> >());

	class_< utility::options::OptionCollection >("OptionCollection",no_init)
			.def("load", &utility::options::OptionCollection::loadVector)
		;

	// class_< core::conformation::ResidueSet, conformation::ResidueSetOP, conformation::ResidueSetCOP >
	//     ( "ResidueSet", init< core::chemical::AtomSetCOP, core::chemical::MMAtomSetCOP >
	//                           (( arg("atom_types_in"), arg("mm_atom_types_in") )) )
	//     .def("read_list_of_residues", &conformation::ResidueSet::read_list_of_residues)
	//   ;

	def("cpp_init", &cpp_init );
	def("cpp_run_tests", & cpp_run_tests ) ;

	def("make_res_set", & make_res_set ) ;
	def("res_set_op_to_ref", &res_set_op_to_ref, return_value_policy< reference_existing_object >() );
	def("res_set_ref_to_op", &res_set_ref_to_op );

	def( "simple_benchmark_smart_ptr", & simple_benchmark_smart_ptr );
	def( "simple_benchmark_ptr",       & simple_benchmark_ptr );
	def( "simple_benchmark_ref",       & simple_benchmark_ref );
	def( "simple_benchmark_by_value",  & simple_benchmark_by_value );


	// class_< core::chemical::AtomSet, core::chemical::AtomSetOP, core::chemical::AtomSetCOP >
	//       ( "AtomSet", init< >() )
	//        .def("read_file", &core::chemical::AtomSet::read_file)
	//   ;
	// class_< core::chemical::MMAtomSet, core::chemical::MMAtomSetOP,core::chemical::MMAtomSetCOP >
	//     ( "MMAtomSet", init< >() )
	//       .def("read_file", &core::chemical::MMAtomSet::read_file)
	//   ;
	// class_< core::scoring::Etable, core::scoring::EtableOP >
	//     ( "Etable", init< core::chemical::AtomSetCOP, core::scoring::EtableOptions const & >() )
	//   ;
	// class_< core::scoring::EtableOptions >( "EtableOptions", init< >() );
	// class_< core::scoring::methods::EnergyMethodManager >
	//     ( "EnergyMethodManager", no_init )
	//   // .def("get_instance", &core::scoring::methods::EnergyMethodManager::get_instance,
	//   //             return_value_policy<reference_existing_object>() )
	//   .def("add_etable",   &core::scoring::methods::EnergyMethodManager::add_etable )
	// ;

	def("EnergyMethodManager_get_instance", &EnergyMethodManager::get_instance,
			return_value_policy<reference_existing_object>() );

	implicitly_convertible<ResidueSetOP,ResidueSetCOP>();
	implicitly_convertible<   AtomSetOP,   AtomSetCOP>();
	implicitly_convertible< MMAtomSetOP, MMAtomSetCOP>();

	enum_< scoring::ScoreType>("ScoreType")
			.value("fa_atr", scoring::fa_atr)
			.value("fa_rep", scoring::fa_rep)
			.value("fa_sol", scoring::fa_sol)
			.value("n_ci_2b_score_types", scoring::n_ci_2b_score_types)
			.value("fa_pair", scoring::fa_pair)
			.value("fa_plane", scoring::fa_plane)
			.value("hbond_sr_bb", scoring::hbond_sr_bb)
			.value("hbond_lr_bb", scoring::hbond_lr_bb)
			.value("hbond_bb_sc", scoring::hbond_bb_sc)
			.value("hbond_sc", scoring::hbond_sc)
			.value("rama", scoring::rama)
			.value("fa_dun", scoring::fa_dun)
			.value("total_score", scoring::total_score)
			.value("n_score_types", scoring::n_score_types)
		.export_values()
		;

	// scope().attr("PHI")   = core::kinematics::PHI;
	// scope().attr("THETA") = core::kinematics::THETA;
	// scope().attr("D")     = core::kinematics::D;
	scope().attr("RAMA")  = scoring::rama;

	class_<Pose>("Pose");
	class_<MoveMap>("MoveMap", init<>() )
		.def( "set_bb", (void ( MoveMap::* )( int const, bool const ) )( &MoveMap::set_bb ), ( arg("seqpos"), arg("setting") ) )
		.def( "set_bb", (void ( MoveMap::* )( bool const ) )( &MoveMap::set_bb ), ( arg("setting") ) )
		.def( "set_chi", (void ( MoveMap::* )( int const,bool const ) )( &MoveMap::set_chi ), ( arg("seqpos"), arg("setting") ) )
		.def( "set_chi", (void ( MoveMap::* )( bool const ) )( &MoveMap::set_chi ), ( arg("setting") ) )
		;
	def("pose_from_pdb",&pose_from_pdb);
	class_<MinimizerOptions>("MinimizerOptions", init< string const &, bool const, bool const, bool const >() );
	class_<AtomTreeMinimizer>("AtomTreeMinimizer", init<>() )
			.def("run", &AtomTreeMinimizer::run )
		;
	class_<ScoreFunction>("ScoreFunction")
			.def("set_weight", &ScoreFunction::set_weight )
		;


}
#endif


int main( int argc, char * argv [] )
{
	//using namespace core;
	std::cerr << "IGNORING argv/argc: " << argv << " " << argc << std::endl;

	cpp_init();

	conformation::ResidueSet & residue_set( *make_res_set() );
	//#ifdef SHEFFLERDEBUG
	std::cerr << "test_wrap.cc:1252 (" << ")" << std::endl;
	//#endif

	core::pose::Pose p;
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:1257 (" << ")" << std::endl;
//#endif
	pose_from_pdb(p,residue_set,"test_in.pdb");
//#ifdef SHEFFLERDEBUG
std::cerr << "test_wrap.cc:1261 (" << ")" << std::endl;
//#endif

	simple_benchmark( residue_set );

	T("TEST") << "SDL:KFJSDL:J\n";
	exit(0);

	/// Output demo
	std::vector<int> A;  A.push_back(1);  A.push_back(2);  A.push_back(3);  A.push_back(5);
	utility::vector1<int> B;  B.push_back(10);  B.push_back(20);  B.push_back(30);  B.push_back(45);
	std::map<int, std::string> M;  M[1]="one";  M[2]="two";  M[3]="1+2";

	T("Demo") << "vector:" << A << " vector1:" << B << " map:" << M << "\n";

	T("Error", 10) << "Some error here!!!\n";
	T("core.pose") << "Some core pose message\n";
	T("core") << "Some core message\n";

	Error() << "Some error test...\n";
	Warning() << "Some warning test...\n";

	exit(0);


	fa_scorefxn_test( residue_set );
	//exit(0);

	//simple_dna_test( residue_set );
	//exit(0);

	simple_benchmark( residue_set );
	exit(0);

	//simple_benchmark( residue_set );
	//exit(0);

	//simple_rotamer_test( residue_set );

	//simple_rotamer_test( residue_set );

	rotamer_trials_test( residue_set );

	pack_rotamers_test( residue_set );
	//exit(0);

	//fa_scorefxn_test( residue_set );

	//test_rama( residue_set );

	simple_conformation_test( residue_set );

	atom_tree_torsion_test( residue_set );

	exit(0);


	//simple_hbond_test( atom_set, residue_set );
	//exit(0);

	//test_dunbrack_io( residue_set );


	//simple_copy_test( residue_set );
	//exit(0);

	//test_rama( residue_set, atom_set );

	// some tests

	simple_min_test( residue_set );


	fa_scorefxn_test( residue_set );

}










































