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
#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/LoopClass.hh>
#include <protocols/frags/TorsionFragment.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hackelec/HackElecEnergy.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/visualize.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/util.hh>
//#include <core/options/after_opts.hh>

#include <core/util/prof.hh> // profiling
#include <core/util/basic.hh>
#include <core/util/SequenceMapping.hh>

#include <core/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>

//silly using/typedef


#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;

//numeric::random::RandomGenerator RG(12323); // <- Magic number, do not change it!!!

using namespace core;
using namespace protocols;

using utility::vector1;

using io::pdb::dump_pdb;


util::Tracer tt( "demo.phil.dna_spec_test", util::t_info );


///////////////////////////////////////////////////////////////////////////////
// some silly helper routines:
//
void
read_list_file( std::string const & filename, utility::vector1< std::string > & l )
{
	std::ifstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message( "cant open list file: "+filename );
	}
	std::string line;
	while ( getline( data, line ) ) {
		std::cout << line << std::endl;
		l.push_back( line );
	}
}


///////////////////////////////////////////////////////////////////////////////
// accumulates into energymap
void
retrieve_residue_pair_energies(
	scoring::EnergyGraph const & energy_graph,
	Size const pos1,
	Size const pos2,
	scoring::EnergyMap & emap
)
{

	using namespace scoring;

	assert( pos2 > pos1 );
	for ( graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node( pos1 )->const_upper_edge_list_begin(),
					irue = energy_graph.get_node( pos1 )->const_upper_edge_list_end();
				iru != irue; ++iru ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
		if ( Size( edge->get_second_node_ind() ) == pos2 ) {
			emap += edge->energy_map();
			break;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
show_clashes(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	Real const clash_threshold
)
{
	using namespace scoring;
	using namespace chemical;

	scorefxn( pose );

	// cached energies object
	Energies & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph & energy_graph( energies.energy_graph() );

	for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		conformation::Residue const & resl( pose.residue( i ) );
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			conformation::Residue const & resu( pose.residue( j ) );

			// the pair energies cached in the link
			EnergyMap const & emap( edge->energy_map());
			Real const clash( emap[ fa_rep ] );
			if ( clash > clash_threshold ) {
				std::cout << "clash: " << resl.name1() << I(4,i) << ' ' << resu.name1() << I(4,j) << F(9,3,clash ) <<
					std::endl;
			}
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
detect_interface_residues(
	pose::Pose const & pose,
	utility::vector1< Size > const & pos_list,
	Real const contact_threshold,
	utility::vector1< bool > & interface
)
{
	Size const nres( pose.total_residue() );

	interface.resize( nres, false );

	for ( Size j=1; j<= nres; ++j ) {
		interface[j] = false;

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( !rsd2.is_protein() ) continue;

		for ( Size n=1; n<= pos_list.size(); ++n ) {
			conformation::Residue const & rsd1( pose.residue( pos_list[n] ) );
			assert( rsd1.is_DNA() );

			if ( rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) ) <= contact_threshold ) {
				interface[ j ] = true;
				break;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
detect_allatom_interface_residues(
	pose::Pose const & pose,
	utility::vector1< Size > const & pos_list,
	Real const contact_threshold,
	utility::vector1< bool > & interface
)
{
	Size const nres( pose.total_residue() );

	interface.resize( nres, false );

	for ( Size j=1; j<= nres; ++j ) {
		interface[j] = false;

		conformation::Residue const & rsd2( pose.residue( j ) );
		if ( !rsd2.is_protein() ) continue;

		for ( Size n=1; n<= pos_list.size(); ++n ) {
			conformation::Residue const & rsd1( pose.residue( pos_list[n] ) );
			assert( rsd1.is_DNA() );

			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {

					if ( rsd1.xyz(ii).distance( rsd2.xyz(jj) ) <= contact_threshold ) {
						interface[ j ] = true;
						break;
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
show_pairing_info(
									pose::Pose const & pose,
									Size const seqpos
									)
{
	using namespace ObjexxFCL::fmt;

	using namespace scoring;
	using namespace scoring::dna;

	DNA_BasePotential const & potential( ScoringManager::get_instance()->get_DNA_BasePotential() );

	BasePartner const & base_partner( retrieve_base_partner_from_pose( pose ) );


	Size const seqpos_partner( base_partner[seqpos] );


	show_base_pair_params( pose.residue(seqpos), pose.residue(seqpos_partner) );

	show_base_step_params( pose.residue(seqpos-1        ), pose.residue(seqpos  ) );
	show_base_step_params( pose.residue(seqpos          ), pose.residue(seqpos+1) );
	show_base_step_params( pose.residue(seqpos_partner-1), pose.residue(seqpos_partner  ) );
	show_base_step_params( pose.residue(seqpos_partner  ), pose.residue(seqpos_partner+1) );

	std::cout << "base_pair_energies: " <<
		F(9,3,potential.base_pair_score( pose.residue( seqpos ), pose.residue( seqpos_partner ) ) ) << std::endl;

	std::cout << "base_step_energies: " <<
		F(9,3, potential.base_step_score( pose.residue( seqpos-1         ), pose.residue(seqpos  ) ) ) <<
		F(9,3, potential.base_step_score( pose.residue( seqpos           ), pose.residue(seqpos+1) ) ) <<
		F(9,3, potential.base_step_score( pose.residue( seqpos_partner-1 ), pose.residue(seqpos_partner  ) ) ) <<
		F(9,3, potential.base_step_score( pose.residue( seqpos_partner   ), pose.residue(seqpos_partner+1) ) ) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
setup_atom_number(
	pose::Pose const & pose,
	id::AtomID_Map< int > & atom_number
)
{
	id::initialize( atom_number, pose );

	int number(0);
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		for ( Size j=1; j<= Size(pose.residue(i).natoms()); ++j ) {
			atom_number[ id::AtomID( j, i ) ] = ++number;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
show_rasmol_hbonds(
	scoring::hbonds::HBondSet const & hbond_set,
	id::AtomID_Map< int > const & atom_number,
	std::ostream & out
)
{
	using namespace scoring::hbonds;

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		if ( hbond_set.allow_hbond(i) ) {
			HBond const & hb( hbond_set.hbond(i) );
			id::AtomID const hatm( hb.don_hatm(), hb.don_res() );
			id::AtomID const aatm( hb.acc_atm() , hb.acc_res() );
			out << "monitor " << atom_number[ hatm ] << ' ' << atom_number[ aatm ] <<
				'\n';
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
wriggle_test(
						 pose::Pose const & start_pose,
						 Size const seqpos
						 )
{
	using namespace conformation;
	using namespace pose;
	using namespace id;

	Real const max_rot( 5.0 ); // degrees

	Residue const &
		prev_rsd( start_pose.residue(seqpos-1) ),
				 rsd( start_pose.residue(seqpos  ) ),
		next_rsd( start_pose.residue(seqpos+1) );

	assert( rsd.is_DNA() && prev_rsd.is_DNA() );

	utility::vector1< Vector > bonds, atoms;

	atoms.push_back( prev_rsd.xyz("O3*") );
	atoms.push_back(      rsd.xyz(  "P") );
	atoms.push_back(      rsd.xyz("O5*") );
	atoms.push_back(      rsd.xyz("C5*") );
	atoms.push_back(      rsd.xyz("C4*") );

	Vector const r( rsd.xyz("O3*") );

	Vector const n1( ( next_rsd.xyz(  "P") -      rsd.xyz("O3*" ) ).normalized() );
	Vector const n2( ( next_rsd.xyz("O5*") - next_rsd.xyz(  "P" ) ).normalized() );

	Real l1( 0.0 );
	utility::vector1< Real > v1(4), v2(4);

	for ( int i=1; i<= 4; ++i ) {
		Vector const bi( ( atoms[i+1] - atoms[i] ).normalized() );
		Vector const ci( atoms[i+1] );
		v1[i] = n1.dot( bi.cross( r - ci ) );
		v2[i] = n2.dot( bi.cross( r - ci ) );
		l1 += v1[i] * v1[i];
	}

	// want to choose e[i] so that e dot v1 = 0 and e dot v2 = 0
	// eg choose two indices randomly
	// choose random values for those indices
	// now choose the unique values of the other two e[i]'s that will give e.v1= e.v2=0
	//
	// even simpler, orthonormalize v1 and v2, then we start with e' and subtract off dot(e',v1)*v1 and
	// dot(e',v2)*v2 and we are golden.
	//

	// normalize v1, compute dot product with v2
	l1 = std::sqrt( l1 );
	Real v1_dot_v2(0.0);
	for ( int i=1; i<= 4; ++i ) {
		v1[i] /= l1;
		v1_dot_v2 += v1[i] * v2[i];
	}

	// now subtract
	Real l2(0.0);
	for ( int i=1; i<= 4; ++i ) {
		v2[i] -= v1_dot_v2 * v1[i];
		l2 += v2[i] * v2[i];
	}
	// normalize v2
	l2 = std::sqrt(l2);
	for ( int i=1; i<= 4; ++i ) v2[i] /= l2;

	// confirm orthogonality
	Real tmp1(0.0), tmp2(0.0), tmp12(0.0);
	for ( int i=1; i<= 4; ++i ) {
		tmp1  += v1[i] * v1[i];
		tmp2  += v2[i] * v2[i];
		tmp12 += v1[i] * v2[i];
	}
	assert( std::abs(tmp1-1)<1e-3 && std::abs(tmp2-1)<1e-3 && std::abs(tmp12)<1e-3 );


	// now try a bunch of moves:
	for ( int r=1; r<= 25; ++r ) {

		for ( int rr=0; rr< 4; ++rr ) {
			bool const subdot1( rr%2 == 1 );
			bool const subdot2( rr/2 == 1 );

		// choose random starting points
		utility::vector1< Real > e(4);
		Real dot1(0.0), dot2(0.0);
		for ( int i=1; i<= 4; ++i ) {
			e[i] = max_rot - 2*numeric::random::uniform()*max_rot;
			dot1 += e[i] * v1[i];
			dot2 += e[i] * v2[i];
		}

		for ( int i=1; i<= 4; ++i ) {
			if ( subdot1 ) e[i] -= dot1 * v1[i];
			if ( subdot2 ) e[i] -= dot2 * v2[i];
		}

		{ //debug
			dot1 = 0; dot2 = 0;
			for ( int i=1; i<= 4; ++i ) {
				dot1 += e[i] * v1[i];
				dot2 += e[i] * v2[i];
			}
			if ( subdot1 ) assert( std::abs( dot1 ) < 1e-3 );
			if ( subdot2 ) assert( std::abs( dot2 ) < 1e-3 );
		}

		// now make torsion angle changes:
		Pose pose;
		pose = start_pose;

		Size const nbb( rsd.n_mainchain_atoms() );
		std::cout << "e[i] ";
		for ( int i=1; i<= 4; ++i ) {

			TorsionID const id( ( i==1 ) ? ( TorsionID( seqpos-1, BB, nbb ) ) : ( TorsionID( seqpos, BB, i-1 ) ) );
			pose.set_torsion( id, pose.torsion(id) + e[i] );
			std::cout << F(9,3,e[i]);
		}
		std::cout << std::endl;

		if ( false ) { // chi
			TorsionID const id( seqpos, CHI, 1 );
			pose.set_torsion( id, pose.torsion(id) + max_rot - 2*numeric::random::uniform()*max_rot );
		}

		std::string tag;
		if ( subdot1 ) tag += "dot1";
		else tag += "____";
		if ( subdot2 ) tag += "dot2";
		else tag += "____";

		dump_pdb( pose, "test_"+tag+"_"+lead_zero_string_of(seqpos,4)+"_"+lead_zero_string_of(r,4)+".pdb" );

		// check dr:
		Vector const r0( pose.residue(seqpos).xyz("O4*") - start_pose.residue(seqpos).xyz("O4*") );
		Vector const r ( pose.residue(seqpos).xyz("O3*") - start_pose.residue(seqpos).xyz("O3*") );

		std::cout << "r: " << tag << F(9,3,r.length()) << F(9,3,r.dot( n1 )) << F(9,3,r.dot( n2 ) ) << std::endl;
		} // rr
	} // r
}



///////////////////////////////////////////////////////////////////////////////
void
dump_hbond_pdb(
	pose::Pose & pose,
	std::string const & filename
)
{
	using namespace optimization;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;

	//
	scoring::hbonds::HBondSet hbond_set;
	pose.update_residue_neighbors();
	scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

	// setup the global atom numbering that would be used for pdb output
	id::AtomID_Map< int > atom_number;
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
show_residue_residue_clashes(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	scoring::ScoreFunction const & scorefxn,
	Real const threshold,
	std::string const & output_tag
	)
{

	using namespace scoring;
	using namespace scoring::methods;
	using namespace scoring::etable;
	using namespace scoring::etable::count_pair;
	using namespace conformation;
	//using namespace scoring::hbonds;

	EtableEnergy const etable_energy
		( *ScoringManager::get_instance()->etable( scorefxn.energy_method_options().etable_type() ),
			scorefxn.energy_method_options() );

	CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

	// assuming only a single bond right now
	// also assuming crossover of 4, should be closest (?) to classic rosetta

	for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
		for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {

			Real weight;
			if ( cpfxn->count( ii, jj, weight ) ) {
				EnergyMap emap;
				Real dsq;
				etable_energy.atom_pair_energy( rsd1.atom(ii), rsd2.atom(jj), weight, emap, dsq );
				if ( emap[ fa_rep ] > threshold ) {
					std::cout << "CLASH: " << F(9,3,emap[ fa_rep ] ) << F(9,3,dsq) <<
						rsd1.name3() << ' ' << rsd1.atom_name(ii) << ' ' <<
						rsd2.name3() << ' ' << rsd2.atom_name(jj) << ' ' <<
						I(4,rsd1.seqpos()) << I(4,rsd2.seqpos() ) << ' ' << output_tag << std::endl;
				}
			}
		}
	}

}


///////////////////////////////////////////////////////////////////////////////
void
show_protein_DNA_interactions(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn
)
{
	using namespace scoring;
	using namespace scoring::methods;
	using namespace scoring::etable;
	using namespace scoring::hackelec;
	using namespace conformation;
	using namespace scoring::hbonds;

	//
	HBondSet hbond_set;
	pose.update_residue_neighbors();
	fill_hbond_set( pose, false, hbond_set );

	hbond_set.sort_by_weighted_energy();

	EtableEnergy const etable_energy
		( *ScoringManager::get_instance()->etable( scorefxn.energy_method_options().etable_type() ),
			scorefxn.energy_method_options() );
	HackElecEnergy const elec_energy( scorefxn.energy_method_options() );

	for ( int i=1; i<= hbond_set.nhbonds(); ++i ) {
		HBond const & hb( hbond_set.hbond(i) );
		Residue const & don_rsd( pose.residue( hb.don_res() ) );
		Residue const & acc_rsd( pose.residue( hb.acc_res() ) );
		if ( don_rsd.is_DNA() || acc_rsd.is_DNA() ) {
// 		if ( ( don_rsd.is_DNA() && acc_rsd.is_protein() ) ||
// 				 ( acc_rsd.is_DNA() && don_rsd.is_protein() ) ) {
			// calculate also the atr,rep,sol and hack_elec energies of interaction, maybe also the sigmoidal electrostatic energy

			Size const dhatm( hb.don_hatm() ), datm( don_rsd.atom_base( dhatm ) ), aatm( hb.acc_atm() );

			// atr,rep,sol: Donor, Hydrogen, Acceptor
			EnergyMap emap;
			Real dsq;
			etable_energy.atom_pair_energy( don_rsd.atom(  datm ), acc_rsd.atom( aatm ), 1.0, emap, dsq );
			etable_energy.atom_pair_energy( don_rsd.atom( dhatm ), acc_rsd.atom( aatm ), 1.0, emap, dsq );
			Real elecE(0.0);
			elecE += elec_energy.eval_atom_atom_hack_elecE( don_rsd.xyz(  datm ), don_rsd.atomic_charge(  datm ),
																											acc_rsd.xyz(  aatm ), acc_rsd.atomic_charge(  aatm ) );
			elecE += elec_energy.eval_atom_atom_hack_elecE( don_rsd.xyz( dhatm ), don_rsd.atomic_charge( dhatm ),
																											acc_rsd.xyz(  aatm ), acc_rsd.atomic_charge(  aatm ) );

			std::cout << "HBOND: " << F(9,3,hb.energy() * hb.weight() ) << ' ' <<
				don_rsd.name1() << I(4,don_rsd.seqpos()) << ' ' << acc_rsd.name1() << I(4,acc_rsd.seqpos()) << ' ' <<
				don_rsd.atom_name( hb.don_hatm() ) << ' ' << acc_rsd.atom_name( hb.acc_atm() ) <<
				F(9,3,hb.energy() ) << F(9,3,hb.weight()) <<
				" energies: " << F(9,3,emap[fa_atr]) << F(9,3,emap[fa_rep]) << F(9,3,emap[fa_sol]) << F(9,3,elecE) << std::endl;
		}
	}

	dump_hbond_pdb( pose, "test_hbond.pdb");

}


////

void
find_dna_rotamers(
 std::ofstream & pdb_lib,
 std::ofstream & dihedral_lib
 )
{
	using utility::vector1;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace conformation;
	using namespace optimization;
	using namespace pose;

	using namespace numeric; // Matrix

	vector1< std::string > const files( options::start_files() );

	Size rot_couter( 0 );

	for ( Size nn=1; nn<= files.size(); ++nn ) {
		Pose pose;
		io::pdb::pose_from_pdb( pose, files[nn] );
		if ( pose.empty() ) continue;

		set_base_partner( pose ); // fills base partner info
		// BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

		for ( Size i=1; i<= pose.total_residue()-1; ++i ) {

			Residue const & rsd     ( pose.residue(i) );

			if ( !rsd.is_DNA() ) continue;
			if ( rsd.is_terminus() ) continue;

			Residue const & prev_rsd( pose.residue(i-1) );
			Residue const & next_rsd( pose.residue(i+1) );

			if ( prev_rsd.is_upper_terminus() || next_rsd.is_lower_terminus() ) continue;

			Pose mini_pose;
			mini_pose.append_residue_by_bond( prev_rsd );
			mini_pose.append_residue_by_bond(      rsd );
 			mini_pose.append_residue_by_bond( next_rsd );

			// kinematics::Stub const takeoff(      rsd.xyz("P"), prev_rsd.xyz("O3*"), prev_rsd.xyz("C3*") );
 			kinematics::Stub const takeoff(      rsd.xyz("C4*"), rsd.xyz("O4*"), rsd.xyz("C1*") );
 			// kinematics::Stub const landing( next_rsd.xyz("P"),      rsd.xyz("O3*"),      rsd.xyz("C3*") );

			id::AtomID_Mask mask;
			id::initialize( mask, mini_pose );

			for ( Size j=1; j<= rsd.natoms(); ++j ) {
				id::AtomID const atom_id(  );
				mask[ id::AtomID( j,2 ) ] = ( ( j<=rsd.last_backbone_atom() )
																			|| j==rsd.first_sidechain_atom() );
#ifndef NDEBUG
				// release build does not compile...
				mini_pose.set_xyz(  id::AtomID( j,2 ), takeoff.global2local( rsd.xyz( j ) ) );
#endif
			}

			for ( Size j=1; j<= next_rsd.natoms(); ++j ) {
#ifndef NDEBUG
				mask[ id::AtomID( j,3 ) ] = j==1;
				mini_pose.set_xyz(  id::AtomID( j,3 ), takeoff.global2local( next_rsd.xyz( j ) ) );
#endif
			}

			pdb_lib << "MODEL     " <<  I(4,++rot_couter) << "\n";
			dump_pdb(mini_pose, pdb_lib, mask);
			pdb_lib << "ENDMDL\n\n";

			dihedral_lib << files[nn]   << ' ' << prev_rsd.name1() << ' '
									 << rsd.name1() << ' ' << next_rsd.name1() << ' '
									 << I( 4,rsd.seqpos() ) << ' ' << rsd.chain();
			dihedral_lib << F(10,3,rsd.mainchain_torsion(6) );
			for (Size i=1; i <= 5; ++i)
				dihedral_lib << F(10,3,rsd.mainchain_torsion(i) );
			dihedral_lib << F(10,3,rsd.chi(1) ) << "\n";


// 			Real const   alpha( rsd.mainchain_torsion(1) );
// 			Real const    beta( rsd.mainchain_torsion(2) );
// 			Real const   gamma( rsd.mainchain_torsion(3) );
// 			Real const   delta( rsd.mainchain_torsion(4) );
// 			Real const epsilon( rsd.mainchain_torsion(5) );
// 			Real const    zeta( rsd.mainchain_torsion(5) );
// 			Real const     chi( rsd.chi(1) );

// 			Vector const o3( rsd.xyz( "O3*" ) ); // look at RotamerSet_.cc dna rotamer code for examples of getting positions




// 			Matrix const transform( landing.M * takeoff.M.transposed() );

// 			Vector const translation( landing.v - takeoff.v );

// 			Vector const local_coordinate_translation( takeoff.global_to_local( landing.v ) );
// 			Vector const local_coordinate_translation( takeoff.global_to_local( next_rsd.xyz("P") ) );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
get_base_pucker(
	conformation::Residue const & rsd,
	std::pair< std::string, std::string > & pucker
)
{

	utility::vector1< std::string > names;
	names.push_back( "C1*" );
	names.push_back( "C2*" );
	names.push_back( "C3*" );
	names.push_back( "C4*" );
	names.push_back( "O4*" );

	utility::vector1< Vector > atoms;
	for ( int i=1; i<= 5; ++i ) {
		atoms.push_back( rsd.xyz( names[i] ) );
	}

	Real mindot = 1000.0;
	for ( int ii=1; ii<= 5; ++ii ) {

		Vector n12 = (( atoms[2]-atoms[1] ).cross( atoms[3]-atoms[2] ) ).normalized();
		Real dot = std::abs( n12.dot( ( atoms[4]-atoms[3] ).normalized() ) );
		if ( dot < mindot ) {
			// get pucker
			Real pucker_dot = n12.dot( ( atoms[5] - Real(0.5) * ( atoms[4] + atoms[1] ) ).normalized() );

			mindot = dot;
			pucker.first = names[5];
			pucker.second = ( pucker_dot < 0.0 ? "endo" : "exxo" );
		}

		atoms.push_back( atoms[1] );
		atoms.erase( atoms.begin() );

		names.push_back( names[1] );
		names.erase( names.begin() );

	}

}


///////////////////////////////////////////////////////////////////////////////
/**
	 Collect stats on intra-dna interactions for the purposes of parametrizing DNA reference energies

	 fa_atr, rep, sol
	 hbond
	 hack_elec

**/

void
intra_dna_stats()
{
	using utility::vector1;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace conformation;
	using namespace optimization;
	using namespace pose;

	ScoreFunction scorefxn;
	scorefxn.energy_method_options().exclude_DNA_DNA( false );


	scorefxn.set_weight( fa_atr, 1.0 );
	scorefxn.set_weight( fa_rep, 1.0 );
	scorefxn.set_weight( fa_sol, 1.0 );
	scorefxn.set_weight( hbond_lr_bb, 1.0 );
	scorefxn.set_weight( hbond_sr_bb, 1.0 );
	scorefxn.set_weight( hbond_bb_sc, 1.0 );
	scorefxn.set_weight( hbond_sc, 1.0 );
	scorefxn.set_weight( hack_elec, 1.0 );


	vector1< std::string > files( options::start_files() );

	for ( Size nn=1; nn<= files.size(); ++nn ) {
		Pose pose;
		io::pdb::pose_from_pdb( pose, files[nn] );
		if ( pose.empty() ) continue;

		set_base_partner( pose ); // fills base partner info
		BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

		scorefxn(pose);

		EnergyGraph const & energy_graph( pose.energies().energy_graph() );

		std::string const tag( lead_zero_string_of( nn, 4 )+".pdb" );

		io::pdb::dump_pdb( pose, tag );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			if ( !rsd.is_DNA() ) continue;

			// base pair energies
			if ( partner[ i ] > i ) {
				EnergyMap emap;
				retrieve_residue_pair_energies( energy_graph, i, partner[i], emap );
				std::cout << "BasePairEnergies: " << rsd.name1() << pose.residue(partner[i]).name1() << ' ' <<
					F(9,3,emap[fa_atr]) << ' ' << F(9,3,emap[fa_rep]) << ' ' << F(9,3,emap[fa_sol]) <<
					' ' << F(9,3,emap[hbond_sc]) << ' ' << F(9,3,emap[ hack_elec ]) << std::endl;
			}

			// base step energies
			if ( partner[i] > i && partner[i+1] != i && partner[i+1] == partner[i]-1 ) {
				EnergyMap emap;
				retrieve_residue_pair_energies( energy_graph, i, i+1, emap );
				retrieve_residue_pair_energies( energy_graph, i, partner[i+1], emap );
				retrieve_residue_pair_energies( energy_graph, i+1, partner[i], emap );
				retrieve_residue_pair_energies( energy_graph, partner[i+1], partner[i], emap );
				std::cout << "BaseStepEnergies: " << rsd.name1() << pose.residue(i+1).name1() << ' ' <<
					F(9,3,emap[fa_atr]) << ' ' << F(9,3,emap[fa_rep]) << ' ' << F(9,3,emap[fa_sol]) <<
					' ' << F(9,3,emap[hbond_sc]) << ' ' << F(9,3,emap[ hack_elec ]) << std::endl;

				std::string dtag1, dtag2;
				{ // i->i+1 dihedrals
					Size const seqpos( i );
					Residue const & rsd1( pose.residue( seqpos   ) );
					Residue const & rsd2( pose.residue( seqpos+1 ) );

					EnergyMap E;
					retrieve_residue_pair_energies( energy_graph, seqpos, seqpos+1, E );

					std::pair< std::string, std::string > p1, p2;
					get_base_pucker( rsd1, p1 );
					get_base_pucker( rsd2, p2 );

					std::ostringstream os;
					os << I(4,seqpos) << " to " << I(4,seqpos+1) <<
						F(9,3,rsd1.chi(1) ) <<
						F(9,3,rsd1.mainchain_torsion(4)) <<
						F(9,3,rsd1.mainchain_torsion(5)) <<
						F(9,3,rsd1.mainchain_torsion(6)) <<
						F(9,3,rsd2.mainchain_torsion(1)) <<
						F(9,3,rsd2.mainchain_torsion(2)) <<
						F(9,3,rsd2.mainchain_torsion(3)) <<
						F(9,3,rsd2.mainchain_torsion(4)) <<
						F(9,3,rsd2.chi(1) ) <<
						F(9,3,E[fa_atr]) << F(9,3,E[fa_rep]) << ' ' <<
						p1.first << '_' << p1.second << ' ' <<
						p2.first << '_' << p2.second;

					std::cout << "dihedrals " << os.str() << '\n';
					dtag1 = os.str();
				}

				{ // partner side i->i+1 dihedrals
					Size const seqpos( partner[i+1] );
					Residue const & rsd1( pose.residue( seqpos   ) );
					Residue const & rsd2( pose.residue( seqpos+1 ) );

					EnergyMap E;
					retrieve_residue_pair_energies( energy_graph, seqpos, seqpos+1, E );

					std::pair< std::string, std::string > p1, p2;
					get_base_pucker( rsd1, p1 );
					get_base_pucker( rsd2, p2 );

					std::ostringstream os;
					os << I(4,seqpos) << " to " << I(4,seqpos+1) <<
						F(9,3,rsd1.chi(1) ) <<
						F(9,3,rsd1.mainchain_torsion(4)) <<
						F(9,3,rsd1.mainchain_torsion(5)) <<
						F(9,3,rsd1.mainchain_torsion(6)) <<
						F(9,3,rsd2.mainchain_torsion(1)) <<
						F(9,3,rsd2.mainchain_torsion(2)) <<
						F(9,3,rsd2.mainchain_torsion(3)) <<
						F(9,3,rsd2.mainchain_torsion(4)) <<
						F(9,3,rsd2.chi(1) ) <<
						F(9,3,E[fa_atr]) << F(9,3,E[fa_rep]) << ' ' <<
						p1.first << '_' << p1.second << ' ' <<
						p2.first << '_' << p2.second;

					std::cout << "dihedrals " << os.str() << '\n';
					dtag2 = os.str();
				}

				show_residue_residue_clashes( rsd, pose.residue(i+1), scorefxn, 0.1, tag+" "+dtag1 );
				show_residue_residue_clashes( rsd, pose.residue(partner[i+1]), scorefxn, 0.1, tag  );
				show_residue_residue_clashes( pose.residue(i+1), pose.residue(partner[i]), scorefxn, 0.1, tag  );
				show_residue_residue_clashes( pose.residue(partner[i+1]), pose.residue(partner[i]), scorefxn, 0.1, tag+" "+dtag2  );
			}
		} // i
	} // files[nn]

}



///////////////////////////////////////////////////////////////////////////////

void
kono_sarai_stats()
{
	using utility::vector1;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace pose;

	vector1< std::string > files( options::start_files() );

	for ( Size nn=1; nn<= files.size(); ++nn ) {
		Pose pose;
		io::pdb::pose_from_pdb( pose, files[nn] );
		if ( pose.empty() ) continue;

// 		set_base_partner( pose ); // fills base partner info
// 		BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

// 		std::string const tag( lead_zero_string_of( nn, 4 )+".pdb" );

// 		io::pdb::dump_pdb( pose, tag );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {

			Residue const & rsd( pose.residue(i) );
			if ( !rsd.is_DNA() ) continue;

			// define coordinate frame
			bool const AG( rsd.aa() == na_ade || rsd.aa() == na_gua );

			Vector const origin( AG ? rsd.xyz("N9") : rsd.xyz("N1") );
			Vector p( AG ? rsd.xyz("C4") : rsd.xyz("C2") );

			Vector const x( ( p - origin ).normalized() );
			Vector const z( get_z_axis( rsd, x ) );
			Vector const y( z.cross( x ) );


			for ( Size j=1; j<= pose.total_residue(); ++j ) {

				Residue const & rsd2( pose.residue(j) );
				if ( !rsd2.is_protein() ) continue;

				Vector const calpha( rsd2.xyz("CA") - origin );
				Real const xx( dot(calpha,x) );
				Real const yy( dot(calpha,y) );
				Real const zz( dot(calpha,z) );

				if ( ( xx >= -13.5 && xx <= 13.5  ) &&
						 ( yy >= -13.5 && yy <= 13.5  ) &&
						 ( zz >=  -6.0 && zz <=  6.0  ) ) {

					std::cout << "CONTACT " << F(9,3,xx) << F(9,3,yy) << F(9,3,zz) <<
						' ' << rsd.name1() << rsd2.name1() << std::endl;
				}
			}
		}
	}
}




///////////////////////////////////////////////////////////////////////////////

void
spec_test(
					pose::Pose const & start_pose,
					scoring::ScoreFunction const & scorefxn,
					utility::vector1< Size > const & pos_list,
					std::string const & output_tag
					)
{
	using namespace pose;
	using namespace devel::dna;
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace optimization;

	Size const nloop( 20 );
	Size const nres( start_pose.total_residue() );



	///////////////////////////////////////////////////////////////////////
	{ // test dna packing

		bool const design( false );


		// setup the task
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( start_pose ));
		kinematics::MoveMap mm;
		{
			//utility::vector1< bool > interface( nres, false );
			//detect_interface_residues( start_pose, pos_list, 12.0, interface );

			task->initialize_from_command_line();

			for ( Size ii = 1; ii <= nres; ++ii ) {
				if ( start_pose.residue(ii).is_protein() ) {
					if ( false ) { //interface[ii] ) {
						task->nonconst_residue_task( ii ).restrict_to_repacking();
						mm.set_chi( ii, true );
					} else {
						task->nonconst_residue_task( ii ).prevent_repacking();
					}
				} else {
					if ( design ) {
						task->nonconst_residue_task( ii ).allow_aa( na_ade );
						task->nonconst_residue_task( ii ).allow_aa( na_thy );
						task->nonconst_residue_task( ii ).allow_aa( na_gua );
						task->nonconst_residue_task( ii ).allow_aa( na_cyt );
					} else {
						task->nonconst_residue_task( ii ).restrict_to_repacking();
					}
				}
			}
		}

		{ // setup residue couplings
			BasePartner const & partner( retrieve_base_partner_from_pose( start_pose ) );
			using namespace pack::rotamer_set;
			RotamerCouplingsOP couplings( new RotamerCouplings() );
			couplings->resize( nres );
			for ( Size i=1;i<= nres; ++i ){
				if ( partner[i] ) {
					(*couplings)[i].first = partner[i];
					(*couplings)[i].second = new conformation::WatsonCrickResidueMatcher();
				}
			}
			task->rotamer_couplings( couplings );
		}

		Pose pose;
		pose = start_pose;

		ScoreFunction scorefxn;
		if ( false ) {
			scorefxn.set_weight( dna_bp, 1.0 );
			scorefxn.set_weight( dna_bs, 1.0 );
		} else {
			scorefxn.set_weight( fa_atr, 0.80 );
			scorefxn.set_weight( fa_rep, 0.44 );
			scorefxn.set_weight( fa_sol, 0.65 );
			scorefxn.set_weight( hbond_lr_bb, 2.0 );
			scorefxn.set_weight( hbond_sr_bb, 2.0 );
			scorefxn.set_weight( hbond_bb_sc, 2.0 );
			scorefxn.set_weight( hbond_sc, 2.0 );
		}

		scorefxn( pose );
		std::cout << "native energies: " << pose.energies().total_energies().show_nonzero() << std::endl;
		show_protein_DNA_interactions( pose, scorefxn );
		show_clashes( pose, scorefxn, 5.0 );
		dump_pdb( pose, "prepack.pdb" );

		for ( Size n=1; n<= pos_list.size(); ++n ) {
			Size const seqpos( pos_list[n] );
			std::cout << "nat pairings: " << seqpos << std::endl;
			show_pairing_info( pose, seqpos );
		}

		scorefxn(pose);
		utility::vector1< std::pair< Real, std::string > > results;
		pack::pack_rotamers_loop( pose, scorefxn, task, nloop, results );

		std::cout << "nat_packed energies: " << pose.energies().total_energies().show_nonzero() << std::endl;
		show_protein_DNA_interactions( pose, scorefxn );
		show_clashes( pose, scorefxn, 5.0 );
		dump_pdb( pose, "postpack.pdb" );
		exit(0);
	}






	/**
		 Go through the positions in pos_list one at a time, mutating to the different bases
		 and showing the energies by component.
	**/

	// setup the task
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( start_pose ));
	kinematics::MoveMap mm;
	{
		utility::vector1< bool > interface( nres, false );
		detect_interface_residues( start_pose, pos_list, 12.0, interface );

		task->initialize_from_command_line();

		for ( Size ii = 1; ii <= nres; ++ii ) {
			if ( start_pose.residue(ii).is_protein() ) {
				if ( interface[ii] ) {
					task->nonconst_residue_task( ii ).restrict_to_repacking();
					mm.set_chi( ii, true );
				} else {
					task->nonconst_residue_task( ii ).prevent_repacking();
				}
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
	}

	{ // analyze native
		Pose pose;
		pose = start_pose;

		scorefxn( pose );
		std::cout << "native energies: " << pose.energies().total_energies().show_nonzero() << std::endl;
		show_protein_DNA_interactions( pose, scorefxn );

		for ( Size n=1; n<= pos_list.size(); ++n ) {
			Size const seqpos( pos_list[n] );
			std::cout << "nat pairings: " << seqpos << std::endl;
			show_pairing_info( pose, seqpos );
		}

		AtomTreeMinimizer().run( pose, mm, scorefxn, MinimizerOptions("dfpmin",0.001,true ) );

		std::cout << "nat_chi-min energies: " << pose.energies().total_energies().show_nonzero() << std::endl;
		show_protein_DNA_interactions( pose, scorefxn );

		pose = start_pose;
		scorefxn(pose);
		utility::vector1< std::pair< Real, std::string > > results;
		pack::pack_rotamers_loop( pose, scorefxn, task, nloop, results );

		std::cout << "nat_packed energies: " << pose.energies().total_energies().show_nonzero() << std::endl;
		show_protein_DNA_interactions( pose, scorefxn );

		AtomTreeMinimizer().run( pose, mm, scorefxn, MinimizerOptions("dfpmin",0.001,true ) );

		std::cout << "nat_packed-chimin energies: " << pose.energies().total_energies().show_nonzero() << std::endl;
		show_protein_DNA_interactions( pose, scorefxn );
	}


	util::prof_reset();
	//BasePartner const & base_partner( retrieve_base_partner_from_pose( start_pose ) );

	for ( Size n=1; n<= pos_list.size(); ++n ) {
		Size const seqpos( pos_list[n] );
		//Size const seqpos_partner( base_partner[seqpos] );

		AA const nat_na( start_pose.residue(seqpos).aa() );
		EnergyMap nat_emap;

		for ( Size r=1;r<=2; ++r ) {

		for ( int nn=first_DNA_aa; nn<= last_DNA_aa; ++nn ) {

			AA na; na = AA(nn);

			if ( ( r == 1 && na != nat_na ) ||
					 ( r == 2 && na == nat_na ) ) continue;

			Pose pose;
			pose = start_pose;

			// make the mutation
			make_base_pair_mutation( pose, seqpos, na );

			// repack at nbring positions
			scorefxn(pose); // hack -- fills 10A nbr graph
			utility::vector1< std::pair< Real, std::string > > results;
			pack::pack_rotamers_loop( pose, scorefxn, task, nloop, results );

			Energy pack_score = scorefxn( pose );

			std::string seq;
			for ( Size i=1; i<= pos_list.size(); ++i ) {
				Size const pos( pos_list[i] );
				if ( pos == seqpos ) {
					seq += uppercased( start_pose.residue(pos).name1() );
				} else {
					seq += pose.residue(pos).name1();
				}
			}

			if ( r == 1 ) {
				nat_emap = pose.energies().total_energies();

				std::cout << "nat_mut pairings: "<< std::endl;

				show_pairing_info( pose, seqpos );
				show_protein_DNA_interactions( pose, scorefxn );


			} else {
				EnergyMap emap( pose.energies().total_energies() );
				emap -= nat_emap;

				std::cout << "mutate: " << output_tag << I(4,seqpos) << ' ' << seq <<
					" from " << start_pose.residue(seqpos).aa() <<
					" to " << na << ' ' << F(9,3,pack_score) << ' ' << emap.show_nonzero() << std::endl;

				show_pairing_info( pose, seqpos );
			}
		}
		} // repeat twice
		util::prof_show();
	}
}


///////////////////////////////////////////////////////////////////////////////
void
idealize_tf_pose( pose::Pose & pose )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace conformation;
	using namespace optimization;
	using namespace id;

	Real const contact_distance( 6.0 );
	Real const atom_pair_sdev( 0.25 );
	Real const coord_sdev( 0.1 );
	Size const window(2); // 5rsds

	Size const nres( pose.total_residue() );

	Size const dna_anchor( 91 );
	{ // setup fold_tree
		assert( nres == 107 );
		kinematics::FoldTree f( nres );
		f.new_jump( 47,  91, 85 );
		f.new_jump( 91, 103, 96 );
		f.reorder( 91 );
		pose.fold_tree(f);
	}

	//dump_pdb( pose, "before.pdb" );
	//insert_ideal_mainchain_bonds( 47, pose.conformation() );
	//dump_pdb( pose, "after.pdb" );
	//exit(0);

	// setup the initial constraints set:


	ConstraintSetOP cst_set( new ConstraintSet() );

	{ // test atompair constraints

		// look for all protein-DNA contacts:
		for ( Size i=1; i<= nres; ++i ) {
			Residue const & i_rsd( pose.residue(i) );
			if ( !i_rsd.is_DNA() ) continue;

			for ( Size j=1; j<= nres; ++j ) {
				Residue const & j_rsd( pose.residue(j) );
				if ( !j_rsd.is_protein() ) continue;

				for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
					for ( Size jj = 1; jj<= j_rsd.natoms(); ++jj ) {

						Real const d( i_rsd.xyz(ii).distance( j_rsd.xyz(jj) ) );
						if ( d < contact_distance ) {
							//std::cout << "add_constraint: " <<
							//	i << ' ' << i_rsd.name() << ' ' << i_rsd.atom_name(ii) << ' ' <<
							//	j << ' ' << j_rsd.name() << ' ' << j_rsd.atom_name(jj) << ' ' << d << ' ' << sdev << std::endl;

							cst_set->add_constraint( new AtomPairConstraint( AtomID(ii,i), AtomID(jj,j),
																															 new HarmonicFunc( d, atom_pair_sdev ) ) );
						}
					} // jj
				} // ii
			} // j,  protein position
		} // i,  DNA position


		for ( Size i=1; i<= nres; ++i ) {
			Residue const & i_rsd( pose.residue(i) );
			if ( !i_rsd.is_protein() ) continue;
			for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
				cst_set->add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,dna_anchor), i_rsd.xyz(ii),
																													 new HarmonicFunc( 0.0, coord_sdev ) ) );
			}
		}
		pose.constraint_set( cst_set );

		ScoreFunction scorefxn;
		//scorefxn.set_weight( atom_pair_constraint, 1.0 );
		scorefxn.set_weight( coordinate_constraint, 1.0 );

		std::cout << "start_score: " << scorefxn( pose ) << std::endl;
		pose.dump_pdb( "start.pdb");

		pose::Pose start_pose;
		start_pose = pose;

		for ( Size n=1; n<= Size(option[ out::nstruct ]); ++n ) { // try idealizing
			pose = start_pose;

			// setup the options
			MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/,
																false /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );
			kinematics::MoveMap mm;

			utility::vector1< Size > pos_list;
			for (Size i=1; i<= nres; ++i ) {
				if ( pose.residue(i).is_protein() ) pos_list.push_back( i );
			}

			int counter(0);
			while ( !pos_list.empty() ) {
				Size const seqpos( pos_list[ static_cast< Size >( numeric::random::uniform() * pos_list.size() + 1 ) ] );
				pos_list.erase( std::find( pos_list.begin(), pos_list.end(), seqpos ) );

				insert_ideal_mainchain_bonds( seqpos, pose.conformation() );

				std::cout << "premin:  " << I(4,seqpos) << F(9,3,CA_rmsd(pose,start_pose)) << F(12,3,scorefxn( pose ) ) <<
					std::endl;

				mm.clear();
				mm.set_chi(seqpos,true);
				for ( Size i=seqpos-window; i<= seqpos+window; ++i ) {
					if ( i>= 1 && i<= nres && pose.residue(i).is_protein() ) {
						mm.set_bb(i,true );
						// disallow proline PHI
						if ( pose.residue(i).aa() == chemical::aa_pro ) mm.set( TorsionID( phi_torsion, BB, i ), false );
					}
				}

				AtomTreeMinimizer().run( pose, mm, scorefxn, options );

				std::cout << "postmin: " << I(4,seqpos) << F(9,3,CA_rmsd(pose,start_pose)) << F(12,3,scorefxn( pose ) ) <<
					std::endl;
				++counter;
				//dump_pdb( pose, "post_min_"+lead_zero_string_of( counter,4 )+"_"+lead_zero_string_of(seqpos,4)+".pdb");
			}
			std::cout << "final: " << I(4,n) << F(9,3,CA_rmsd(pose,start_pose)) << F(12,3,scorefxn( pose ) ) <<
				std::endl;

			dump_pdb( pose, "idl_"+lead_zero_string_of( n, 4 )+".pdb" );
		}
		exit(0);




		pose.set_phi( 47, pose.phi(47) + 15 );

		std::cout << "end_score: " << scorefxn( pose ) << std::endl;
		std::cout << "end_score2: " << scorefxn( pose ) << std::endl;

		pose.dump_pdb( "end.pdb");

		// try minimizing:

		{
			using namespace optimization;

			kinematics::MoveMap mm;
			mm.set_bb( 45, true );
			mm.set_bb( 46, true );
			mm.set_bb( 48, true );
			mm.set_bb( 49, true );

			// setup the options
			MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/,
																true /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );

			AtomTreeMinimizer().run( pose, mm, scorefxn, options );
			std::cout << "min_score: " << scorefxn( pose ) << std::endl;

			pose.dump_pdb( "min.pdb");

		}


	} // scope


	{ // test coordinate constraints on all protein positions


	}



}





///////////////////////////////////////////////////////////////////////////////
void
zif268_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;
	using namespace pose;
	using namespace scoring;

	Size const motif_begin( 90 );
	Size const motif_size( 3 );

	std::string const filename( "input/1aay.pdb" ); //option[ phil::s ]+"_subset.pdb" );

	Pose pose;
	io::pdb::pose_from_pdb( pose, filename );

	io::pdb::dump_pdb( pose, "test_zif.pdb" );

	//assert( pose.total_residue() == 32 );

	scoring::dna::set_base_partner( pose );


	// setup scoring function
	//
	utility::vector1< std::string > weights_files;
	if ( option[ weights_tag ]() != "none" ) {
		weights_files.push_back( option[ weights_tag ] );
	} else {
		read_list_file( option[ weights_tag_list ], weights_files );
	}



	for ( Size j=1; j<= weights_files.size(); ++j ) {

		ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( weights_files[j] ) );

		// dna specific mods
		scorefxn->set_weight( fa_pair, 0.0 );
		scorefxn->set_weight( hack_elec, option[ Whack_elec ] );
		scorefxn->set_weight( dna_bp, option[ Wdna_bp ] );
		scorefxn->set_weight( dna_bs, option[ Wdna_bs ] );

		devel::dna::packing_specificity_test( pose, *scorefxn, motif_begin, motif_size, "dfpmin", 0.001, true,
																					option[ output_tag ]+weights_files[j], !option[ fast ] );

	}


}




///////////////////////////////////////////////////////////////////////////////
void
bzip_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;
	using namespace pose;
	using namespace scoring;

	Size const motif_begin( 21 );
	Size const motif_size( 3 );

	std::string const filename( options::start_file()+"_subset.pdb" );

	Pose pose;
	io::pdb::pose_from_pdb( pose, filename );

	assert( pose.total_residue() == 32 );

	scoring::dna::set_base_partner( pose );


	// setup scoring function
	//
	utility::vector1< std::string > weights_files;
	if ( option[ weights_tag ]() != "none" ) {
		weights_files.push_back( option[ weights_tag ] );
	} else {
		read_list_file( option[ weights_tag_list ], weights_files );
	}



	for ( Size j=1; j<= weights_files.size(); ++j ) {

		ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( weights_files[j] ) );

		// dna specific mods
		scorefxn->set_weight( fa_pair, 0.0 );
		scorefxn->set_weight( hack_elec, option[ Whack_elec ] );
		scorefxn->set_weight( dna_bp, option[ Wdna_bp ] );
		scorefxn->set_weight( dna_bs, option[ Wdna_bs ] );

		devel::dna::packing_specificity_test( pose, *scorefxn, motif_begin, motif_size, "dfpmin", 0.001, true,
																					weights_files[j] );

	}



}

///////////////////////////////////////////////////////////////////////////////
void
endo_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;
	using namespace pose;
	using namespace scoring;

	Size const motif_begin( 7 );
	Size const motif_size( 3 );

	std::string const filename( options::start_file()+"_subset.pdb" ); // -s cmd input

	Pose pose;
	io::pdb::pose_from_pdb( pose, filename );

	assert( pose.total_residue() == 133 );

	io::pdb::dump_pdb( pose, "test.pdb" );

	scoring::dna::set_base_partner( pose ); // finds base pairs


	// setup scoring function
	//
	utility::vector1< std::string > weights_files; // declare weights_files vector
	if ( option[ weights_tag ]() != "none" ) { // read in weights from file
		weights_files.push_back( option[ weights_tag ] );
	} else {
		read_list_file( option[ weights_tag_list ], weights_files );
	}



	for ( Size j=1; j<= weights_files.size(); ++j ) {

		ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( weights_files[j] ) );

		// dna specific mods
		scorefxn->set_weight( fa_pair, 0.0 );
		scorefxn->set_weight( hack_elec, option[ Whack_elec ] );
		scorefxn->set_weight( dna_bp, option[ Wdna_bp ] );
		scorefxn->set_weight( dna_bs, option[ Wdna_bs ] );

		devel::dna::packing_specificity_test( pose, *scorefxn, motif_begin, motif_size, "dfpmin", 0.001, true,
																					weights_files[j] );

		utility::vector1< int > motif_positions;
		for ( int i=1; i<= 14; ++i ) {
			if ( i!=10 ) motif_positions.push_back(i);
		}

		devel::dna::packing_specificity_test_fast( pose, *scorefxn, motif_positions, option[ nloop ],
																							 option[ min_type ], option[ minimize_tolerance ],
																							 option[ post_minimize ], weights_files[j] /* tag */ );
	}
}



///////////////////////////////////////////////////////////////////////////////
void
cst_test()
{
	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/1aay.pdb");
	idealize_tf_pose( pose );
}

///////////////////////////////////////////////////////////////////////////////
void
tf_specificity_test(
	std::string const & tf_name,
	std::string wts_tag,
	std::string const & output_tag,
	Real const Whack_elec,
	Real const Wdna_bp,
	Real const Wdna_bs
)
{
	using namespace core::options;
	using namespace core::options::OptionKeys::dna::specificity;
	using namespace scoring;
	using namespace scoring::dna;
	using namespace pose;
	using namespace conformation;

	//std::string const home_directory( "/Users/pbradley/" );
	std::string const home_directory( "/home/pbradley/" );

	// read the motifs file
	// get the pdbfile name, and the motif positions to be repacked
	//
	Real const frequency_threshold( 0.75 );
	Real const contact_dis2( 4.5 * 4.5 );
	//Real const contact_dis2( -999 );

	std::string pdb_file;
	Size dna_chain;
	utility::vector1< Size > motif_positions;
	utility::vector1< char > motif_pdb_seq; // for debugging
	{
		std::string const filename( home_directory+"/dat/dna/am_data/motif_pdb_matches/"+tf_name+".motif" );
		std::ifstream data( filename.c_str() );
		assert( data.good() );
		std::string line,tag,pdb;
		getline( data, line );

		{ // parse header
			std::istringstream l(line);
			l >> tag >> tag >> pdb >> tag >> tag >> dna_chain;
			assert( !l.fail() );
			pdb_file = home_directory+"/dat/dna/am_data/pdb_original/"+pdb+".pdb";
		}
		while ( getline(data,line) ) {
			std::istringstream l(line);
			char seq;
			int pos;
			Real highest_frequency, crystal_frequency;
			l >> tag >> tag >> pos >> seq >> crystal_frequency >> highest_frequency;
			assert( !l.fail() );
			if ( highest_frequency >= frequency_threshold ) {
				assert( std::abs( crystal_frequency - highest_frequency )<1e-3 ); // this will fail, tmp hack
				std::cout << "motif position: " << line << std::endl;
				motif_positions.push_back( pos );
				motif_pdb_seq.push_back( seq );
			}
		}
		data.close();
	}


	// read the pdb file
	Pose pose;
	utility::vector1< Size > new_motif_positions; // in pose's numbering system
	{
		pose::Pose pdb_pose;
		io::pdb::pose_from_pdb( pdb_pose, pdb_file );

		if ( pdb_pose.total_residue() < 1 ) {
			std::cout << "pdb io failed: " << pdb_file << std::endl;
			return;
		}

		for ( Size i=1; i<= motif_positions.size(); ++i ) {
			assert( pdb_pose.chain_sequence(dna_chain)[ motif_positions[i]-1 ] == motif_pdb_seq[i] );
		}

		set_base_partner( pdb_pose );
		BasePartner const & partner( retrieve_base_partner_from_pose( pdb_pose ) );

		// find DNA partner for this dna_chain and the protein chain(s) that interact with it...
		std::set< Size > chains;
		chains.insert( dna_chain );

		for ( Size k=1; k<= motif_positions.size(); ++k ) {
			Size const i( pdb_pose.conformation().chain_begin( dna_chain ) - 1 + motif_positions[k] );
			assert( pdb_pose.sequence()[ i-1 ] == motif_pdb_seq[k] );

			//
			if ( partner[i] ) chains.insert( pdb_pose.chain( partner[i] ) );

			// check for protein contacts
			Residue const & i_rsd( pdb_pose.residue(i) );
			for ( Size j=1; j<= pdb_pose.total_residue(); ++j ) {
				Residue const & j_rsd( pdb_pose.residue(j) );
				if ( !j_rsd.is_protein() ) continue;
				bool contact( false );
				for ( Size ii = 1; ii<= i_rsd.natoms() && !contact; ++ii ) {
					for ( Size jj = 1; jj<= j_rsd.natoms() && !contact; ++jj ) {
						if ( i_rsd.xyz(ii).distance_squared( j_rsd.xyz(jj) ) < contact_dis2 ) contact = true;
					} // jj
				} // ii
				if ( contact ) chains.insert( pdb_pose.chain( j ) );
			} // j protein position

		} // i (k) motif position


		// now add the chains to the new pose:
		for ( std::set< Size >::const_iterator it = chains.begin(); it != chains.end(); ++it ) {
			Size const c( *it );
			Size const chain_begin( pdb_pose.conformation().chain_begin(c) );

			if ( c == dna_chain ) {
				int motif_offset( pose.total_residue() );
				for ( Size k=1; k<= motif_positions.size(); ++k ) {
					new_motif_positions.push_back( motif_positions[k] + motif_offset );
				}
			}

			for ( Size i= chain_begin, i_end= pdb_pose.conformation().chain_end(c); i<= i_end; ++i ) {
				if ( i == chain_begin && !pose.empty() ) {
					pose.append_residue_by_jump( pdb_pose.residue(i), 1 );
					pose.conformation().insert_chain_ending( pose.total_residue()-1 );
				} else {
					pose.append_residue_by_bond( pdb_pose.residue(i) );
				}
			}

		}
	} // scope

	set_base_partner( pose );



	dump_pdb( pose, "test.pdb" );


	// standard packer wts
	if ( wts_tag == "none" ) {
		if ( option[ soft_rep ]  ) {
			wts_tag = SOFT_REP_WTS;
		} else {
			wts_tag = STANDARD_WTS;
		}
	}


	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( wts_tag ) );

	// dna specific mods
	scorefxn->set_weight( fa_pair, 0.0 );
	scorefxn->set_weight( hack_elec, Whack_elec );
	scorefxn->set_weight( dna_bp, Wdna_bp );
	scorefxn->set_weight( dna_bs, Wdna_bs );

	//(*scorefxn)(pose);
	//scorefxn->accumulate_residue_total_energies( pose );
	//pose.energies().show( std::cout );

	{ // make a debugging pdb
		std::string filename( tf_name+"_motif.pdb" );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) { // only write the file once
			data.close();
			std::ofstream out( filename.c_str() );
			out << "load pdb inline\nselect " << new_motif_positions[1];
			for ( Size i=2; i<= new_motif_positions.size(); ++i ) {
				out << ", " << new_motif_positions[i];
			}
			out << "\ndefine motif selected\nwireframe 75\nrestrict not hydrogen\nexit\n";
			dump_pdb( pose, out );
		} else {
			data.close();
		}
	}

	// HACK
	if ( false ) {
		kinematics::FoldTree f( pose.fold_tree() );
		f.new_jump( 209, 211, 210 );
		pose.fold_tree(f);
		wriggle_test( pose, 210 );
		return;
	}


	spec_test( pose, *scorefxn, new_motif_positions, tf_name+" "+output_tag );
	return;


	devel::dna::packing_specificity_test_fast( pose, *scorefxn, new_motif_positions, option[ nloop ],
																						 option[ min_type ], option[ minimize_tolerance ],
																						 option[ post_minimize ], output_tag );
}

void
tf_specificity_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;

	utility::vector1< std::string > tfs( options::start_files()), weights_files;

	if ( option[ weights_tag ]() != "none" ) {
		weights_files.push_back( option[weights_tag] );
	} else {
		read_list_file( option[ weights_tag_list ], weights_files );
	}

	for ( Size i=1; i<= tfs.size(); ++i ) {
		for ( Size j=1; j<= weights_files.size(); ++j ) {
			std::string const outtag( tfs[i] + "_" + weights_files[j] + "_" + option[output_tag] );

			tf_specificity_test( tfs[i], weights_files[j], outtag, option[ Whack_elec ], option[ Wdna_bp ],
													 option[ Wdna_bs ] );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
dna_specificity_test(
	std::string const & filename,
	int const motif_begin,
	int const motif_size,
	std::string wts_tag,
	std::string const & output_tag,
	Real const Whack_elec,
	Real const Wdna_bp,
	Real const Wdna_bs
)
{
	using namespace core::options;
	using namespace core::options::OptionKeys::dna::specificity;
	using namespace scoring;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, filename );

	if ( pose.total_residue() < 1 ) {
		std::cout << "pdb io failed: " << filename << std::endl;
		return;
	}

	//int const motif_begin( 21 );
	//int const motif_size( 3 );

	dna::set_base_partner( pose );

	// standard packer wts
	if ( wts_tag == "none" ) {
		if ( option[ soft_rep ]  ) {
			wts_tag = SOFT_REP_WTS;
		} else {
			wts_tag = STANDARD_WTS;
		}
	}

	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( wts_tag ) );

	// dna specific mods
	scorefxn->set_weight( fa_pair, 0.0 );
	scorefxn->set_weight( hack_elec, Whack_elec );
	scorefxn->set_weight( dna_bp, Wdna_bp );
	scorefxn->set_weight( dna_bs, Wdna_bs );

	if ( option[ fast ] ) {
		devel::dna::packing_specificity_test_fast( pose, *scorefxn, motif_begin, motif_size, option[ nloop ],
																							 option[ min_type ], option[ minimize_tolerance ],
																							 option[ post_minimize ], output_tag );
	} else {
		devel::dna::packing_specificity_test( pose, *scorefxn, motif_begin, motif_size, option[ min_type ],
																					option[ minimize_tolerance ], option[ post_minimize ], "" );
	}

// 	{ // hacking
// 		ScoreFunction scorefxn;
// 		scorefxn.set_weight( fa_atr, 0.8 );
// 		scorefxn.set_weight( fa_rep, 0.4 );
// 		scorefxn.set_weight( fa_sol, 0.6 );
// 		scorefxn.set_weight( fa_intra_rep, 0.00004 );

// 		devel::dna::packing_specificity_test_fast( pose, scorefxn, motif_begin, motif_size, option[ min_type ],
// 																							 option[ post_minimize ] );
// 	}

// 	{ // hacking
// 		ScoreFunction scorefxn;
// 		scorefxn.set_weight( fa_atr, 0.8 );
// 		scorefxn.set_weight( fa_rep, 0.4 );
// 		scorefxn.set_weight( fa_sol, 0.6 );
// 		scorefxn.set_weight( fa_intra_rep, 0.004 );

// 		devel::dna::packing_specificity_test_fast( pose, scorefxn, motif_begin, motif_size, option[ min_type ],
// 																							 option[ post_minimize ] );
// 	}

// 	{ // hacking
// 		ScoreFunction scorefxn;
// 		scorefxn.set_weight( fa_atr, 0.8 );
// 		scorefxn.set_weight( fa_rep, 0.4 );
// 		scorefxn.set_weight( fa_sol, 0.6 );
// 		scorefxn.set_weight( hbond_lr_bb, 1.2 );
// 		scorefxn.set_weight( hbond_sr_bb, 1.2 );
// 		scorefxn.set_weight( hbond_sc, 1.2 );
// 		scorefxn.set_weight( hbond_bb_sc, 1.4 );
// 		devel::dna::packing_specificity_test_fast( pose, scorefxn, motif_begin, motif_size, option[ min_type ],
// 																							 option[ post_minimize ] );
// 	}

// 	{ // hacking
// 		ScoreFunction scorefxn;
// 		scorefxn.set_weight( fa_atr, 0.8 );
// 		scorefxn.set_weight( fa_rep, 0.4 );
// 		scorefxn.set_weight( fa_sol, 0.6 );

// 		devel::dna::packing_specificity_test_fast( pose, scorefxn, motif_begin, motif_size, option[ min_type ],
// 																							 option[ post_minimize ] );
// 	}


// 	{ // hacking
// 		ScoreFunction scorefxn;
// 		scorefxn.set_weight( fa_atr, 0.8 );
// 		scorefxn.set_weight( fa_rep, 0.4 );
// 		scorefxn.set_weight( fa_sol, 0.6 );
// 		scorefxn.set_weight( fa_dun, 0.5 );

// 		devel::dna::packing_specificity_test_fast( pose, scorefxn, motif_begin, motif_size, option[ min_type ],
// 																							 option[ post_minimize ] );
// 	}

// 	{ // hacking
// 		ScoreFunction scorefxn;
// 		scorefxn.set_weight( fa_atr, 0.8 );
// 		scorefxn.set_weight( fa_rep, 0.4 );
// 		scorefxn.set_weight( fa_sol, 0.6 );
// 		scorefxn.set_weight( p_aa_pp, 0.5 );
// 		scorefxn.set_weight( ref, 0.25 );

// 		devel::dna::packing_specificity_test_fast( pose, scorefxn, motif_begin, motif_size, option[ min_type ],
// 																							 option[ post_minimize ] );
// 	}
}

///////////////////////////////////////////////////////////////////////////////
void
loop_modeling_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace id;
	using namespace io::pdb;
	using namespace pose;
	using namespace protocols::loops;
	using namespace protocols::frags;


	// read the pose
	Pose pose;
	io::pdb::centroid_pose_from_pdb( pose, options::start_file() );

	dump_pdb( pose, "premap.pdb" );

	{ // hacking

		for ( int na=first_DNA_aa; na <= last_DNA_aa; ++na ) {
			Real max_d( 0.0 );

			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd( pose.residue(i) );
				if ( rsd.aa() != AA(na) ) continue;

				Vector const & nbr_xyz( rsd.nbr_atom_xyz() );

				for ( Size j=1; j<= rsd.natoms(); ++j ) {
					Real const d( nbr_xyz.distance( rsd.xyz(j) ) );
					if ( d > max_d ) max_d = d;
				}
			}
			std::cout << "MAX_D: " << AA(na) << ' ' << max_d << std::endl;
		}
		//exit(0);
	}

	// read the alignment file
	util::SequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( option[ phil::align_file ], source_seq, target_seq, mapping );


	// may not represent the entire pose source sequence
	std::string const pose_seq( pose.sequence() );

	if ( source_seq != pose_seq ) {
		/// source sequence from align file does not cover entire pdb file sequence

		if ( pose_seq.find( source_seq ) == std::string::npos ) {
			utility_exit_with_message( "alignfile source sequence not contained in pose sequence" );
		}
		Size const offset( pose_seq.find( source_seq ) );

		std::string source_nterm_seq, target_nterm_seq;
		for ( Size i=1; i<= offset; ++i ) {
			// this is a residue in the input pdb that's not accounted for in the alignment
			// if its protein, unalign it
			// if its dna, align it
			Residue const & rsd( pose.residue( i ) );

			mapping.insert_residue( i );
			source_nterm_seq += rsd.name1();

			if ( rsd.is_DNA() ) {
				target_nterm_seq += rsd.name1();
				Size const target_pos( target_nterm_seq.size() );
				mapping.insert_target_residue( target_pos );
				mapping[i] = target_pos;
			} else {
				assert( rsd.is_protein() );
				// ligand etc could be tricky when reading an alignment from a file
			}
		}

		source_seq = source_nterm_seq + source_seq;
		target_seq = target_nterm_seq + target_seq;

		assert( pose_seq.find( source_seq ) == 0 );

		while( source_seq.size() < pose_seq.size() ) {
			assert( mapping.size1() == source_seq.size() );

			Size const pos( source_seq.size() + 1 );
			Residue const & rsd( pose.residue( pos ) );

			char const n1( rsd.name1() );
			source_seq += n1;
			mapping.push_back( 0 );

			if ( rsd.is_DNA() ) {
				target_seq += n1;
				mapping.insert_target_residue( target_seq.size() );
				mapping[pos] = target_seq.size();
			}
		}

		assert( pose_seq == source_seq );
	}

	//mapping.size2( target_seq.size() );

	std::cout << pose.sequence() << '\n' << target_seq << '\n';

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		std::cout << "M1: " << I(4,i) << ' ' << pose.residue(i).name() << I(4,mapping[i]);
		if ( mapping[i] ) {
			std::cout << ' ' << target_seq[ mapping[i]-1 ];
		}
		std::cout << std::endl;
	}

	protocols::loops::apply_sequence_mapping( pose, target_seq, mapping );


	dump_pdb( pose, "start.pdb" );



	Loops loops;

	{ // figure out where the loops are
		util::SequenceMapping m( mapping );
		m.reverse();

		// debug:
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			std::cout << "M2: " << I(4,i) << ' ' << pose.residue(i).name() << I(4,m[i]) << std::endl;
		}

		bool in_loop( false );
		Size loop_begin( 0 );

		for ( Size i=1; i<= m.size1(); ++i ) {
			if ( m[i] == 0 ) {
				if ( !in_loop ) {
					loop_begin = i;
					in_loop = true;
				}
			} else {
				if ( in_loop ) {
					loops.add_loop( loop_begin, i-1, i-1, 0.0, 0, true );
					in_loop = false;
				}
			}
		}

		//loops.add_loop(  82,  89,  84, 0, true );
		//loops.add_loop( 100, 107, 103, 0, true );
	}

	/**
		 fragment setup:
		 1. pick 3mer library from vall
		 2. derive 1mer libary from 3mers


	**/
	utility::vector1<int> frag_sizes; //( option[ OptionKeys::loops::frag_sizes ] );



	frag_sizes.push_back(3);
	frag_sizes.push_back(1);

	//FileVectorOption frag_files( option[OptionKeys::loops::frag_files] );
	//assert( frag_sizes.size() == frag_files.size() );

	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > frag_libs;
	std::map< Size, bool > frag_libs_init;


	{ // scope
		using namespace protocols::frags;
		{ // 3mers
			Size const frag_size( 3 );
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
			// gets 3mers from vall:
			setup_loops_fragment_libraries( pose.sequence(), loops, frag_libs, frag_libs_init );
		}

		{ // 1mers
			// gets 1mers from 3mers
			Size const frag_size = 1;
			Size const prev_size = 3;
			protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
			frag_libs_init.insert( std::make_pair( frag_size , false ) );
			frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );

			protocols::frags::TorsionFragmentLibraryOP    prev_lib_op( frag_libs.find( prev_size)->second );

			frag_libs_init[ frag_size ] = frag_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
		}

		// confirm success
		std::cout << frag_libs_init[1] << frag_libs_init[3] << std::endl;
		assert( frag_libs_init[1] && frag_libs_init[3] );
	}

	// for graphics:
	protocols::viewer::add_conformation_viewer( pose.conformation(), "loops_pose" );

	util::prof_reset();

	perturb_loops_with_ccd( pose, loops, frag_libs );

	util::prof_show();

	{ // switch to fullatom mode, preserve DNA sidechain conformations
		Pose cen_pose;
		cen_pose = pose;

		protocols::loops::switch_to_residue_type_set( pose, FA_STANDARD );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const &     rsd(     pose.residue(i) );
			Residue const & src_rsd( cen_pose.residue(i) );
			if ( rsd.is_DNA() ) {
				assert( rsd.natoms() == src_rsd.natoms() );
				for ( Size j=1; j<= rsd.natoms(); ++j ) {
					if ( rsd.atom_is_backbone(j) ) assert( rsd.xyz(j).distance_squared( src_rsd.xyz( rsd.atom_name(j) ) )<1e-2 );
					pose.set_xyz( AtomID(j,i), src_rsd.xyz( rsd.atom_name(j) ) );
				}
			}
		}
	}


	util::prof_reset();
	refine_loops_with_ccd( pose, loops );

	util::prof_show();

	dump_pdb( pose, "final.pdb" );

}


///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
void
luxr_setup()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
// 	using namespace scoring;
 	using namespace scoring::dna;
 	using namespace conformation;
 	using namespace chemical;
// 	using namespace optimization;
// 	using namespace id;
// 	using namespace io::pdb;
 	using namespace pose;
// 	using namespace protocols::loops;
// 	using namespace protocols::frags;
// 	using namespace protocols::moves;
// 	using namespace scoring::dna; //a

	// read the pose
	Pose start_pose;
	io::pdb::pose_from_pdb( start_pose, options::start_file() );

	Pose pdb_pose;
	pdb_pose = start_pose;


	// read the alignment file
	util::SequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( option[ phil::align_file ], source_seq, target_seq, mapping );

	protocols::loops::extend_sequence_mapping( start_pose /*const*/, mapping, source_seq, target_seq );

	assert( source_seq == start_pose.sequence() && mapping.size1() == source_seq.size() &&
					mapping.size2()==target_seq.size() );

	protocols::loops::apply_sequence_mapping( start_pose, target_seq, mapping );

	assert( start_pose.sequence() == target_seq );

	Size const nres( start_pose.total_residue() ); //a

	// now retrieve the DNA chi angles from pdb pose
	// take coordinates exactly if there's a sequence match -- actually that's already done
	//

	{
		util::SequenceMapping m( mapping );
		Size const nres( start_pose.total_residue() );

		assert( nres == m.size2() );
		m.reverse();
		assert( nres == m.size1() );

		for ( Size i=1; i<= start_pose.total_residue(); ++i ) {
			if ( m[i] && start_pose.residue(i).is_DNA() ) start_pose.set_chi( 1, i, pdb_pose.chi( 1, m[i] ) );
		}
	}

	// a
	int n1,n2,n4;
	// calculate foldtree numbers
	core::Real contact_distance( 6.0 );
	for ( Size i=1; i<= nres; ++i ) {
		Residue const & i_rsd( start_pose.residue(i) );
		if ( !i_rsd.is_DNA() ) continue;

		for ( Size j=1; j<= nres; ++j ) {
			Residue const & j_rsd( start_pose.residue(j) );
			if ( !j_rsd.is_protein() ) continue;

			for ( Size ii = 1; ii<= i_rsd.natoms(); ++ii ) {
				for ( Size jj = 1; jj<= j_rsd.natoms(); ++jj ) {

					core::Real const d( i_rsd.xyz(ii).distance( j_rsd.xyz(jj) ) );
					if ( d < contact_distance ) {
						contact_distance = d;
						n1 = j;
						n2 = i;
					}
				}
			}
		}
	}

	Size const n3( start_pose.conformation().chain_end(1) );
	Size const n5( start_pose.conformation().chain_end(2) );

	{ // helpful error messages
		Size const nchains( start_pose.conformation().num_chains() );
		if ( nchains != 3 ) {
			utility_exit_with_message( "Expected exactly three chains in the pdb file!" );
		}
		if ( !( start_pose.residue(1).is_protein() &&
						start_pose.residue(n3+1).is_DNA() &&
						start_pose.residue(n5+1).is_DNA() ) ) {
			utility_exit_with_message( "Expected three chains: first protein, then dna1, then dna2" );
		}
	}


	// get the base pair of n2
	set_base_partner( start_pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( start_pose ) );
	n4 = partner[n2];

	Size dna_jump(0);
	{
		kinematics::FoldTree f( start_pose.total_residue() );
		f.new_jump( n1, n2, n3 );
		f.new_jump( n2, n4, n5 );
		f.reorder( n2 );
		start_pose.fold_tree(f);
		dna_jump = 1;
		std::cout << start_pose.fold_tree() << std::endl;
	}


	for ( Size i=0; i<= nres; ++i ) {
		if ( i==0 || i == nres || start_pose.chain(i) != start_pose.chain(i+1) ) {
			if ( i>0    ) make_variant_residue( start_pose, UPPER_TERMINUS, i   );
			if ( i<nres ) make_variant_residue( start_pose, LOWER_TERMINUS, i+1 );
		}
	}

	io::pdb::dump_pdb( start_pose, option[ out::file::o ] );
	exit(0);

}

///////////////////////////////////////////////////////////////////////////////
void
calc_protein_DNA_rmsd(
											pose::Pose const & pose,
											pose::Pose const & reference_pose,
											Real & ca_rmsd,
											Real & interface_ca_rmsd,
											Real & interface_allatom_rmsd
											)
{
	Size const nres( pose.total_residue() );

	utility::vector1< int > pos_list;
	utility::vector1< bool > interface;
	for ( Size i=1; i<= nres; ++i ) {
		if ( pose.residue(i).is_DNA()) pos_list.push_back(i);
	}

	detect_allatom_interface_residues( reference_pose, pos_list, 4.5, interface );

	ca_rmsd = 0.0;
	interface_ca_rmsd = 0.0;
	interface_allatom_rmsd = 0.0;

	int protres(0), intres(0), intatoms(0);
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) ), ref_rsd( reference_pose.residue(i) );
		if ( !rsd.is_protein() ) continue;

		Real const ca_dis2( rsd.xyz("CA").distance_squared( ref_rsd.xyz("CA") ) );
		++protres;
		ca_rmsd += ca_dis2;
		if ( interface[i] ) {
			++intres;
			interface_ca_rmsd += ca_dis2;
			for ( Size j=1; j<= rsd.natoms(); ++j ) {
				++intatoms;
				interface_allatom_rmsd += rsd.xyz(j).distance( ref_rsd.xyz( rsd.atom_name(j) ) );
			}
		}
	}

	ca_rmsd = std::sqrt( ca_rmsd/ protres );
	interface_ca_rmsd = std::sqrt( interface_ca_rmsd / intres );
	interface_allatom_rmsd = std::sqrt( interface_allatom_rmsd / intatoms );
}

///////////////////////////////////////////////////////////////////////////////
void
luxr_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;
 	using namespace scoring;
 	using namespace scoring::dna;
 	using namespace pose;

	Size const nstruct( option[ out::nstruct] ); // how many structures to make

	if ( !( option[ motif_begin  ].user() ) || !( option[ motif_size   ].user() ) ||
			 !( option[ moving_jump  ].user() ) ) {
		utility_exit_with_message( "motif_begin and motif_size are required args");
	}
	Size const m_begin( option[ motif_begin ] );
	Size const m_size ( option[ motif_size ] );

	// read the pose
	Pose start_pose;
	io::pdb::pose_from_pdb( start_pose, options::start_file() );

	set_base_partner( start_pose );

	// for graphics: recycle a single pose for simulations
	Pose pose;
	pose = start_pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "rb_min_pose" );

	///////////////////////////
	// setup the score function:
	scoring::ScoreFunctionOP scorefxn;
	if ( option[ score_function ].user() ) {
		scorefxn = new ScoreFunction();
		scorefxn->initialize_from_file( option[ score_function ] );
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( STANDARD_WTS );
		// dna specific mods
		scorefxn->set_weight( fa_pair, 0.0 );
		scorefxn->set_weight( hack_elec, option[ Whack_elec ] ); // was 0.5
		scorefxn->set_weight(    dna_bp, option[ Wdna_bp    ] );
		scorefxn->set_weight(    dna_bs, option[ Wdna_bs    ] );
	}


	// create the protocol object
	devel::dna::RB_Relax rb_relax( *scorefxn, option[ moving_jump ] );

	for ( Size n=1; n<= nstruct; ++n ) {

		// recover the starting conformation
		pose = start_pose;

		// randomize the dna positions
		if ( option[ randomize_motif ] ) devel::dna::randomize_motif_sequence( pose, m_begin, m_size );

		util::prof_reset();

		// run the protocol
		rb_relax.apply( pose );

		util::prof_show();

		std::string const filename( option[ output_tag ] + "final" + lead_zero_string_of( n, 4 ) + ".pdb" );

		Real ca_rmsd, interface_ca_rmsd, interface_allatom_rmsd;

		calc_protein_DNA_rmsd( pose, start_pose, ca_rmsd, interface_ca_rmsd, interface_allatom_rmsd );

		tt << "final_score: " << ca_rmsd << ' ' << interface_ca_rmsd << ' ' << interface_allatom_rmsd << ' ' <<
			pose.sequence().substr(m_begin-1, m_size ) << ' ' <<
			start_pose.sequence().substr(m_begin-1, m_size ) << ' ' <<
			(*scorefxn)(pose) << ' ' << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) <<
			' ' << filename << '\n';

		pose.dump_pdb( filename );

	} // nstruct

}


///////////////////////////////////////////////////////////////////////////////
void
atom_vdw_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace scoring;
	using namespace scoring::constraints;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace id;
	using namespace pose;

	utility::vector1< std::string > pdb_files( options::start_files() );

	AtomTypeSet const & atom_types( *chemical::ChemicalManager::get_instance()->atom_type_set("centroid_dna") );
	Size const natomtypes( atom_types.n_atomtypes() );

	FArray2D< Real > closest( natomtypes, natomtypes, 999 );

	for ( Size n=1; n<= pdb_files.size(); ++n ) {
		Pose pose;
		std::cout << "try to read: " << pdb_files[n] << std::endl;

		io::pdb::centroid_pose_from_pdb( pose, pdb_files[n] );

		Size const nres( pose.total_residue() );

		// look at protein-dna closest-approach distances
		for ( Size i=1; i<= nres; ++i ) {
			Residue const & rsd1( pose.residue(i) );
			if ( !rsd1.is_protein() ) continue;
			for ( Size j=1; j<= nres; ++j ) {
				Residue const & rsd2( pose.residue(j) );
				if ( !rsd2.is_DNA() ) continue;
				for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
					Size const iitype( rsd1.atom_type_index( ii ) );
					for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
						Size const jjtype( rsd2.atom_type_index( jj ) );
						Real const dis2( rsd1.xyz(ii).distance_squared( rsd2.xyz(jj) ) );
						if ( dis2 < closest( iitype, jjtype ) ) {
							closest( iitype, jjtype ) = dis2;
							closest( jjtype, iitype ) = dis2;
						}
					}
				}
			}
		}

		dump_pdb(pose,"test.pdb");
	}

	// dump a new atomvdw file
	AtomVDW const & atom_vdw( ScoringManager::get_instance()->get_AtomVDW() );

	std::ofstream out( "atom_vdw_tmp.txt" );

	for ( Size i=1; i<= natomtypes; ++i ) {
		out << atom_types[i].name();
		utility::vector1< Real > const & i_atom_vdw( atom_vdw(i) );
		for ( Size j=1; j<= natomtypes; ++j ) {
			if ( closest(i,j) < 998 ) {
				out << F(12,6,closest(i,j) );
			} else {
				out << F(12,6,i_atom_vdw[j] );
			}
		}
		out << '\n';
	}

	out.close();

	exit(0);



}

///////////////////////////////////////////////////////////////////////////////
void
dna_specificity_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;

	utility::vector1< std::string > pdb_files( options::start_files() ), weights_files;

	if ( option[ weights_tag ]() != "none" ) {
		weights_files.push_back( option[weights_tag] );
	} else {
		read_list_file( option[ weights_tag_list ], weights_files );
	}

	for ( Size i=1; i<= pdb_files.size(); ++i ) {
		for ( Size j=1; j<= weights_files.size(); ++j ) {
			std::string const outtag( option[output_tag] + "_" + string_of(i) + "_" + string_of(j) );
			std::cout << "FILES: " << outtag << ' ' << pdb_files[i] << ' ' << weights_files[j] << std::endl;

			dna_specificity_test( pdb_files[i], option[ motif_begin ], option[ motif_size ], weights_files[j],
														outtag, option[ Whack_elec ], option[ Wdna_bp ], option[ Wdna_bs ] );
		}
	}


}


///////////////////////////////////////////////////////////////////////////////
void
dump_dna_kinemage()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, options::start_file() );

	kinematics::dump_pose_kinemage( "tmp.kin", pose );
	exit(0);

}


///////////////////////////////////////////////////////////////////////////////
void
read_write()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, options::start_file() );
	pose.dump_pdb( option[ out::file::o ] );

}

///////////////////////////////////////////////////////////////////////////////
void
analyze_interface_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, options::start_file() );

	Real bsasa14, bsasa5;
	utility::vector1< id::AtomID > bdon, bacc;
	devel::dna::analyze_interface_sasa( pose, option[ dna::specificity::moving_jump ], bsasa14, bsasa5, bdon, bacc );

	{ // SILLY
		scoring::ScoreFunction sf;
		sf.set_weight( scoring::hbond_sc, 1.0 );

		sf(pose);
		//dump_hbond_pdb( pose, "test_hbonds.pdb" );

		scoring::hbonds::HBondSet hbond_set;
		pose.update_residue_neighbors();
		scoring::hbonds::fill_hbond_set( pose, false, hbond_set );

		// setup the global atom numbering that would be used for pdb output
		id::AtomID_Map< int > atom_number;
		setup_atom_number( pose, atom_number );

		std::ofstream out( "hbond_unsat.pdb" );
		out << "load pdb inline\n";
		show_rasmol_hbonds( hbond_set, atom_number, out );
		// now show unsat hbond donors and acceptors
		out << "select none\ndefine donors selected\n";
		for ( Size i=1; i<= bdon.size(); ++i ) {
			out << "select donors or atomno=" << atom_number[ bdon[i] ] << "\ndefine donors selected\n";
		}
		out << "select none\ndefine acceptors selected\n";
		for ( Size i=1; i<= bacc.size(); ++i ) {
			out << "select acceptors or atomno=" << atom_number[ bacc[i] ] << "\ndefine acceptors selected\n";
		}

		out << "select donors\ncpk 100\nselect acceptors\ncpk 100\nexit\n";
		dump_pdb( pose, out );
		out.close();
	}
}

///////////////////////////////////////////////////////////////////////////////
void
water_test()
{
	using namespace pose;
	using namespace id;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace core::options;
	using namespace core::options::OptionKeys;


	Pose pdb_pose;
	io::pdb::pose_from_pdb( pdb_pose, options::start_file() );

	// 1st delete waters that arent with 12A of protein AND dna
	Pose new_pose;

	{ // scope
		Pose const pose( pdb_pose );
		//pose = pdb_pose;

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );

			if ( rsd.is_polymer() ) {
				if ( rsd.is_lower_terminus() ) {
					new_pose.append_residue_by_jump( rsd, 1 );
				} else {
					new_pose.append_residue_by_bond( rsd );
				}
			} else {
				assert( rsd.name() == "TP3" );

				// see if we're a nbr water
				bool protein_nbr( false ), dna_nbr( false );
				Vector const & xyz( rsd.nbr_atom_xyz() );

				for ( Size j=1; j<= pose.total_residue(); ++j ) {
					Residue const & rsd2( pose.residue(j) );

					if ( rsd2.is_polymer() ) {
						if ( xyz.distance_squared( rsd2.nbr_atom_xyz() ) < 100 ) {
							if ( rsd2.is_protein() ) {
								protein_nbr = true;
								if ( dna_nbr ) break;
							} else {
								dna_nbr = true;
								if ( protein_nbr ) break;
							}
						}
					}
				}
				if ( protein_nbr && dna_nbr ) {
					new_pose.append_residue_by_jump( rsd, 1 );
				}
			} // water
		} // i

		new_pose.conformation().chains_from_termini();
	} // scope


	// try packing waters
	{
		Pose pose( new_pose );

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		kinematics::MoveMap mm;
		task->initialize_from_command_line();

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).name() == "TP3" ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
				assert( task->pack_residue( ii ) );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
				assert( !task->pack_residue( ii ) );
			}
		}

		ScoreFunction scorefxn;
		scorefxn.set_weight( fa_atr, 0.80 );
		scorefxn.set_weight( fa_rep, 0.44 );
		scorefxn.set_weight( fa_sol, 0.65 );
		scorefxn.set_weight( hbond_lr_bb, 2.0 );
		scorefxn.set_weight( hbond_sr_bb, 2.0 );
		scorefxn.set_weight( hbond_bb_sc, 2.0 );
		scorefxn.set_weight( hbond_sc, 2.0 );

		scorefxn( pose );
		utility::vector1< std::pair< Real, std::string > > results;

		pack::pack_rotamers_loop( pose, scorefxn, task, 25, results );

	} // scope

}



///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	std::string const mode( option[ dna::specificity::mode ].value() );
	if ( mode == "luxr" ) {
		luxr_test();
		exit(0);
	} else if ( mode == "luxr_setup" ) {
		luxr_setup();
		exit(0);
	} else if ( mode == "read_write" ) {
		read_write();
		exit(0);
	} else if ( mode == "intra_dna_stats" ) {
		intra_dna_stats();
		exit(0);
	} else if ( mode == "analyze_interface" ) {
		analyze_interface_test();
		exit(0);
	} else if ( mode == "water" ) {
		water_test();
		exit(0);
	} else if ( mode == "zif268" ) {
		zif268_test();
		exit(0);
	}



// 	// Open output files ....
// 	std::ofstream pdb_lib( "sugar_based_dna_rotamers.pdb" ),
// 		dihedral_lib ( "sugar_based_dna_rotamers.dih" );
// 	if ( !pdb_lib.is_open() || !dihedral_lib.is_open() )
// 		;

	exit(0); // add new mode strings

// 	utility::vector1< conformation::ResidueOP > new_rotamers;
// 	std::ifstream lib_file;
// 	lib_file.open( "/home/cyanover/projects/dna/tmp/VQ_DNA_64_rotamers_7_dof" );
// 	load_dna_rotamers( lib_file, new_rotamers );


// 	find_dna_rotamers(pdb_lib, dihedral_lib);

	tf_specificity_test();
	exit(0);


// 	// Close output files ....
// 	pdb_lib.close();
// 	dihedral_lib.close();

// 	exit(0);

	zif268_test();

	kono_sarai_stats();

	exit(0);

	loop_modeling_test();
	exit(0);

	atom_vdw_test();
	exit(0);

	//dump_dna_kinemage();
	//exit(0);

	{ //
		// Open output files ....
		std::ofstream pdb_lib( "dna_rotamers.pdb" ),
			dihedral_lib ( "dna_rotamers.dih" );
		if ( !pdb_lib.is_open() || !dihedral_lib.is_open() )
			;

		find_dna_rotamers(pdb_lib, dihedral_lib);

		// Close output files ....
		pdb_lib.close();
		dihedral_lib.close();

		exit(0);
	}

	zif268_test();
	exit(0);

	bzip_test();
	exit(0);

	cst_test();
	exit(0);

	endo_test();
	exit(0);

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	// initialize option and random number system
	core::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}

		// 		// now do specificity calculation

		// 		devel::dna::packing_specificity_test_fast( pose, *scorefxn, motif_begin, motif_size, 50, // nloop
		// 																							 "none", 1000.0, //option[ min_type ], option[ minimize_tolerance ],
		// 																							 false, filename+"yesWdna", true, false );

		// 		scorefxn->set_weight( dna_bp, 0.0 );
		// 		scorefxn->set_weight( dna_bs, 0.0 );

		// 		devel::dna::packing_specificity_test_fast( pose, *scorefxn, motif_begin, motif_size, 50, // nloop
		// 																							 "none", 1000.0, //option[ min_type ], option[ minimize_tolerance ],
		// 																							 false, filename+"noWdna", true, false );

		// 		scorefxn->set_weight( dna_bp, option[ Wdna_bp ] );
		// 		scorefxn->set_weight( dna_bs, option[ Wdna_bs ] );

// 	for ( Size i = 1; i <= frag_sizes.size(); ++i ) {
// 		Size const frag_size = Size(frag_sizes[i]);
// 		protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
// 		frag_libs_init.insert( std::make_pair(frag_size, false ) );//frag_lib_op->read_file( frag_files[i], frag_size, 3 ) ) );
// 		frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
// 	}

// 	Size prev_size(10000);
// 	protocols::frags::TorsionFragmentLibraryOP prev_lib_op(0);
// 	for ( std::map<Size, bool>::const_reverse_iterator it = frag_libs_init.rbegin(),
// 					it_end = frag_libs_init.rend(); it != it_end; it++ ) {
// 		Size const frag_size( it->first );
// 		bool const frag_lib_init( it->second );
// 		assert( frag_size < prev_size );
// 		if ( (!frag_lib_init) && prev_lib_op ) {
// 			std::cout << "set up " << frag_size << "-mer library from " << prev_size << "-mer library" << std::endl;
// 			protocols::frags::TorsionFragmentLibraryOP current_lib_op( frag_libs.find(frag_size)->second );
// 			frag_libs_init[frag_size] = current_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
// 		}
// 		prev_size = frag_size;
// 		prev_lib_op = frag_libs[frag_size];
// 		std::cout << "frag_libs_init: " << frag_size << " " << frag_libs_init[frag_size] << std::endl;
// 	}

// 	setup_loops_fragment_libraries( pose.sequence(), loops, frag_libs, frag_libs_init );
