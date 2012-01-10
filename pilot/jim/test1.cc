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
#include <protocols/frags/VallData.hh>
#include <protocols/frags/TorsionFragment.hh>
#include <devel/dna/protocols.hh>

#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/Methods.hh>

#include <protocols/moves/BackboneMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/OutputMovers.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/RigidBodyMover.hh>
//#include <protocols/moves/rigid_body_moves.hh>
#include <protocols/moves/TrialMover.hh>

#include <protocols/loops/ccd_closure.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/viewer/viewers.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>

#include <core/util/prof.hh> // profiling
#include <core/util/CacheableData.hh> // profiling

#include <core/util/SequenceMapping.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
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
#include <core/kinematics/visualize.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/util.hh>//option.hh>
//#include <core/options/after_opts.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>

#include <core/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>

//REMOVE LATER!
//#include <utility/io/izstream.hh>


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

static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

using namespace core;
using namespace protocols;

using utility::vector1;

using io::pdb::dump_pdb;


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
show_intrachain_energies(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn
)
{
	using namespace scoring;
	using namespace chemical;

	scorefxn( pose );

	//Size const nres( pose.total_residue() );
	Size const nchain( pose.conformation().num_chains() );

	FArray2D< EnergyMap > chain_energies( nchain, nchain );

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

			chain_energies( resl.chain(), resu.chain() ) += edge->energy_map();
		}
	}

	for ( Size i=1; i<= nchain; ++i ) {
		for ( Size j=i; j<= nchain; ++j ) {
			std::cout << "chainE:  " << i << ' ' << j;
			for ( Size k=1; k<= n_score_types; ++k ) {
				if ( scorefxn[ ScoreType(k) ] ) {
					std::cout << ' ' << ScoreType(k) << ' ' << F(9,3,chain_energies(i,j)[ ScoreType(k) ]);
				}
			}
			std::cout << std::endl;
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
test_dunbrack_io()
{
	//using namespace core;
	using namespace scoring;
	using namespace scoring::dunbrack;

	RotamerLibrary const & rot_lib
		( ScoringManager::get_instance()->get_RotamerLibrary() );

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		rot_lib.rotamer_energy( pose.residue(i) );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
test_gb()
{
	//using namespace core;
	using namespace scoring;
	using namespace pose;
	using namespace conformation;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, start_file() );

	ScoreFunction scorefxn;
	scorefxn.set_weight( gb_elec, 1.0 );

	scorefxn( pose );


	// show radii
	GenBornPoseInfo const & gb_info( static_cast< GenBornPoseInfo & >( pose.data().get( util::GEN_BORN_POSE_INFO ) ) );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		GenBornResidueInfo const & gb( gb_info.residue_info( i ) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			std::cout << "GBINFO: " << I(4,i) << ' ' << rsd.name3() << ' ' << rsd.atom_name(j) <<
				F(9,2,gb.born_radius(j) ) << F(9,2,gb.atomic_radius(j) ) << std::endl;
		}
	}

	LREnergyContainerCOP lrec = pose.energies().long_range_container( methods::gen_born_lr );
	assert ( !lrec->empty() );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		std::cout << "GBPAIR: " << ii << ' ' << ii  << ' ' <<
			F(9,2,pose.energies().onebody_energies( ii )[ gb_elec ] ) << std::endl;

		for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_upper_neighbor_iterator_begin( ii ),
							rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
			Size const jj = rni->upper_neighbor_id();
			EnergyMap emap;
			rni->retrieve_energy( emap );
			std::cout << "GBPAIR: " << ii << ' ' << jj  << ' ' << F(9,2,emap[gb_elec]) << std::endl;
		}
	}


	dump_pdb( pose, "new_gb.pdb" );

	exit(0);

	/// now pack
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

	task->initialize_from_command_line();
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ii%2 == 0 ) task->nonconst_residue_task( ii ).restrict_to_repacking();
		else task->nonconst_residue_task( ii ).prevent_repacking();
	}

	pack::pack_rotamers( pose, scorefxn, task);

}

///////////////////////////////////////////////////////////////////////////////
void
test_rama()
{
	using namespace scoring;
	std::cout << "test rama" << std::endl;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	ScoreFunction scorefxn;
	scorefxn.set_weight( rama , 0.2 );

	scorefxn( pose );

	std::cout << "rama succeeded" << std::endl;
}




///////////////////////////////////////////////////////////////////////////////
void
test_scorefxn_io()
{
	using namespace scoring;
	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	// standard packer wts
	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( STANDARD_WTS ) );
	scorefxn->set_weight( fa_rep, 0.0 );
	(*scorefxn)(pose);

	std::cout << *scorefxn;

	/// soft rep packer wts
	ScoreFunctionOP scorefxn2( ScoreFunctionFactory::create_score_function( SOFT_REP_WTS ) );
	(*scorefxn2)(pose);

	std::cout << *scorefxn2;

	/// score12 w/ std packer wts
	ScoreFunctionOP score12( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	(*score12)(pose);

	std::cout << *score12;
}

///////////////////////////////////////////////////////////////////////////////
void
simple_rotamer_test()
{
	using namespace conformation;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	ResidueOP ile_rsd( pose.residue(3).clone() );
	ResidueOP pro_rsd( pose.residue(70).clone() );

	assert( pro_rsd->name() == "PRO" && ile_rsd->name() == "ILE" );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		if ( rsd.is_terminus() ) continue;

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
///


///////////////////////////////////////////////////////////////////////////////
///

void
set_fullatom_flag_test()
{
	using namespace io::pdb;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	using namespace protocols::frags;
	using namespace protocols::moves;

	{ // if we dont already have a fullatom pose
		pose::Pose pose;
		centroid_pose_from_pdb( pose, "input/test_in.pdb" );

		// now convert to fullatom:
		ResidueTypeSetCAP fullatom_residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			Residue const & rsd( pose.residue(i) );
			for ( Size j=1; j<=rsd.type().variant_types().size(); j++ ) {
				std::cout << rsd.type().variant_types()[j] << std::endl;
			}
			// get all residue types with same AA
			ResidueTypeCAPs const & rsd_types( fullatom_residue_set->aa_map( rsd.aa() ) );
			ResidueOP new_rsd( 0 );
			// now look for a rsdtype with same variants
			for ( Size j=1; j<= rsd_types.size(); ++j ) {
				ResidueType const & new_rsd_type( *rsd_types[j] );
				if ( rsd.type().variants_match( new_rsd_type ) ) {
					new_rsd = ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
					break;
				}
			}
			assert( new_rsd ); // really should print error msg

			pose.replace_residue( i, *new_rsd, false );
		}

		dump_pdb( pose, "fullatom_test1.pdb" );
	} // scope


	/// other way: if we already have a fullatom_pose
	{

		pose::Pose fullatom_pose;
		pose_from_pdb( fullatom_pose, "input/test_in.pdb" );

		pose::Pose pose;
		centroid_pose_from_pdb( pose, "input/test_in.pdb" );

		// now convert to fullatom:
		// only tricky thing -- cutpoint residue types might not be preserved...
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			pose.replace_residue( i, fullatom_pose.residue(i), true /* orient_backbone */ );
		}

		dump_pdb( pose, "fullatom_test2.pdb" );

	} // scope

}

void
simple_loop_modeling_test()
{
	using namespace io::pdb;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace id;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	using namespace protocols::frags;
	using namespace protocols::moves;

	using core::Real;
	Real const init_phi  ( -150.0 );
	Real const init_psi  (  150.0 );
	Real const init_omega(  180.0 );

	// get rid of this crap when the MonteCarlo changes are reverted
	pose::PoseOP pose_ptr( new pose::Pose() );
	pose::Pose & pose( *pose_ptr );
	centroid_pose_from_pdb( pose, "input/test_in.pdb" );

	std::cout << pose.sequence() << std::endl;

	util::SequenceMapping mapping;
	std::string source_seq, target_seq;
	read_alignment_file( start_file(), source_seq, target_seq, mapping );

	util::SequenceMapping start_mapping( mapping );

	assert( mapping.size1() == source_seq.size() && mapping.size2() == target_seq.size() );

	std::cout << "start mapping: " << std::endl;
	mapping.show();


	// alternate approach:

	// first delete all unaligned residues
	while ( !mapping.all_aligned() ) {
		for ( Size i=1; i<= mapping.size1(); ++i ) {
			if ( !mapping[i] ) {
				pose.conformation().delete_polymer_residue( i );
				mapping.delete_residue( i );
				break; // because the numbering is screwed up now, start the loop again
			}
		}
	}

	// now convert sequence of aligned positions
	ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
	{
		for ( Size i=1; i<= mapping.size1(); ++i ) {
			char const new_seq( target_seq[ mapping[i]-1 ] ); // strings are 0-indexed
			// will fail if .select(...) returns empty list
			ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( new_seq ).exclude_variants().select( rsd_set )[1] );
			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue(i), pose.conformation() ) );
			pose.replace_residue( i, *new_rsd, false );
		}
	}

	// setup fold_tree
	{
		assert( pose.total_residue() == mapping.size1() );
		kinematics::FoldTree f( pose.fold_tree() );
		for ( Size i=1; i< pose.total_residue(); ++i ) {
			assert( mapping[i] );
			if ( mapping[i+1] != mapping[i]+1 ) {
				f.new_jump( i, i+1, i );
			}
		}
		pose.fold_tree(f);
		std::cout << "FOldtree: " << f << std::endl;
	}

	// add terminal residues
	// nterm
	while ( mapping[ 1 ] != 1 ) {
		int const aligned_pos( mapping[1] - 1 );
		char const new_seq( target_seq[ aligned_pos-1 ] ); // 0-indexed
		ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( new_seq ).exclude_variants().select( rsd_set )[1] );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
		pose.conformation().prepend_polymer_residue_before_seqpos( *new_rsd, 1, true );
		pose.set_omega( 1, 180.0 );
		mapping.insert_aligned_residue( 1, aligned_pos );
	}
	// cterm
	while ( mapping[ mapping.size1() ] != mapping.size2() ) {
		int const seqpos( mapping.size1() + 1 );
		int const aligned_pos( mapping[seqpos-1] + 1 );
		char const new_seq( target_seq[ aligned_pos-1 ] ); // 0-indexed
		ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( new_seq ).exclude_variants().select( rsd_set )[1] );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
		pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, seqpos-1, true );
		pose.set_omega( seqpos-1, 180.0 );
		mapping.insert_aligned_residue( seqpos, aligned_pos );
	}

	// now fill in the breaks
	{
		for ( Size cut=1; cut<= Size(pose.fold_tree().num_cutpoint()); ++cut ) {
			while ( mapping[ pose.fold_tree().cutpoint( cut )+1] != Size(pose.fold_tree().cutpoint( cut )+1 )) {
				Size const cutpoint( pose.fold_tree().cutpoint( cut ) );
				assert( mapping[cutpoint] == cutpoint ); // we've fixed everything up til here

				// add at the nterm of the loop
				int const aligned_pos( cutpoint+1 );
				char const new_seq( target_seq[ aligned_pos - 1 ] ); // 0-indexed
				ResidueTypeCAP new_rsd_type( ResidueSelector().set_name1( new_seq ).exclude_variants().select( rsd_set )[1] );
				ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
				pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, cutpoint, true );
				pose.set_omega( cutpoint, 180.0 );
				mapping.insert_aligned_residue( cutpoint + 1, aligned_pos );
			} // while missing residues
			if ( true ) {
				int const cutpoint( pose.fold_tree().cutpoint( cut ) );
				make_variant_residue( pose, CUTPOINT_LOWER, cutpoint   );
				make_variant_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
			}
		} // cut=1,ncutpoints
	} // scope

	assert( mapping.is_identity() );
	assert( pose.sequence() == target_seq );

	dump_pdb( pose, "loops_start.pdb" );

	// setup allow move
	Size const nres( pose.total_residue() );
	MoveMap mm;
	start_mapping.reverse();
	assert( start_mapping.size1() == nres );
	for ( Size i=1; i<= nres; ++i ) {
		if ( ( !start_mapping[i] ) || // unaligned
				 ( i>   1 && start_mapping[i] != start_mapping[i-1]+1 ) ||
				 ( i<nres && start_mapping[i] != start_mapping[i+1]-1 ) ) {
			std::cout << "moving: " << i << ' ' << start_mapping[i] << std::endl;
			mm.set_bb( i, true );
			if ( !start_mapping[i] ) {
				pose.set_phi  ( i, init_phi   );
				pose.set_psi  ( i, init_psi   );
				pose.set_omega( i, init_omega );
			}
		}
	}


	////////////////////////
	// pick fragments
	// read vall:
	VallData vall( option[ phil::vall_file ] );

	TorsionFragmentLibrary lib;
	Size const frag_size(3);

	lib.resize( nres - frag_size + 1 );
	for ( Size i=1; i<= nres-frag_size+1; ++i ) {
		std::string const frag_seq( target_seq.substr(i-1,3) );
		vall.get_frags( 200, frag_seq, "---", 1.0, 0.0, false, false, true, lib[i] );
	}


	// setup scoring function
	// get rid of this crap when MonteCarlo changes are reverted
	core::scoring::ScoreFunctionOP scorefxn_ptr( new ScoreFunction() );
	ScoreFunction & scorefxn( *scorefxn_ptr );

	scorefxn.set_weight( scoring::env, 1.0 );
	scorefxn.set_weight( scoring::pair, 1.0 );
	scorefxn.set_weight( scoring::cbeta, 1.0 );
	scorefxn.set_weight( scoring::vdw, 1.0 );
	scorefxn.set_weight( scoring::chainbreak, 1.0 );

	for ( Size i=1; i<= nres; ++i ) {
		std::cout << "natoms " << i << ' ' << pose.residue(i).natoms() << std::endl;
	}

	scorefxn( pose );
	pose.energies().show( std::cout );

	{ // testing
		pose::Pose new_pose;
		new_pose = pose;
		Energies E( pose.energies() );

		std::ofstream out1( "out1" );
		std::ofstream out2( "out2" );
		std::ofstream out3( "out3" );

		pose.energies().show( out1 );
		new_pose.energies().show( out2 );
		E.show( out3 );

		out1.close();
		out2.close();
		out3.close();

	}

	{ // hacking
		FoldTree f( nres );
		f.new_jump( 1, 30, pose.fold_tree().cutpoint(1) );
		pose.fold_tree(f);
	}


	// fragment moves
	MonteCarlo mc( *pose_ptr, *scorefxn_ptr, 2.0 );

	protocols::viewer::add_conformation_viewer( pose.conformation() );
	protocols::viewer::add_monte_carlo_viewer( mc );

	util::prof_reset();
	for ( int i=1; i<= 10; ++i ) {
		for ( int ii=1; ii<= 100; ++ii ) {
			// choose an insertion position
			int pos;
			while ( true ) {
				pos = static_cast< int > ( ( nres - frag_size + 1 ) * RG.uniform() + 1 );
				bool allowed( true );
				for ( Size k=0; k< frag_size; ++k ) {
					if ( !( mm.get( TorsionID( pos+k, BB,   phi_torsion ) ) &&
									mm.get( TorsionID( pos+k, BB,   psi_torsion ) ) &&
									mm.get( TorsionID( pos+k, BB, omega_torsion ) ) ) ) {
						allowed = false;
						break;
					}
				}
				if ( allowed ) break;
			}

			// choose a fragment
			Size const nfrags( lib[ pos ].size() );
			int const nn( static_cast< int >( nfrags * RG.uniform() + 1 ) );

			// make move
			lib[pos][nn].insert( pose, pos );

			std::cout << "score_before: " << scorefxn( pose ) << std::endl;

			mc.boltzmann( *pose_ptr );

			std::cout << "score_after: " << scorefxn( pose ) << ' ' << mc.last_accepted_score() << ' ' <<
				mc.lowest_score() << std::endl;

		}
		mc.recover_low( *pose_ptr );

		util::prof_show();
	}

	util::prof_show();
}

// for graphics
void*
simple_loop_modeling_test_wrapper( void* )
{
	simple_loop_modeling_test();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
void
dna_deriv_test_old()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace optimization;
	using namespace options;
	using namespace pose;

	using namespace options::OptionKeys;
	using namespace io::pdb;
	using namespace protocols::frags;
	using namespace protocols::moves;
	using namespace scoring::dna;

	using core::Real;

	//
	Pose pdb_pose;
	pose_from_pdb( pdb_pose, start_file()); // eg 1qn4
	Size const nres( pdb_pose.total_residue() );

	set_base_partner( pdb_pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( pdb_pose ) );

	// create a mini-pose with just dna
	Pose pose, jump_pose;

	std::string const jump_anchor_atom( " C2 " );
	for ( Size chain1_begin=1; chain1_begin<= nres; ++chain1_begin ) {
		if ( partner[ chain1_begin ] ) {
			// found first dna residue
			int const chain1( pdb_pose.chain( chain1_begin ) );

			for ( Size j=chain1_begin; ( j<= nres && partner[j] && pdb_pose.chain(j) == chain1 ); ++j ) {
				Residue const &         rsd( pdb_pose.residue(            j ) );
				Residue const & partner_rsd( pdb_pose.residue( partner[ j ] ) );

				Size const rsd_seqpos( j - chain1_begin + 1 );
				if ( rsd_seqpos == 1 ) {
					pose.append_residue_by_bond( rsd );
					jump_pose.append_residue_by_bond( rsd );
				} else {
					pose.append_polymer_residue_after_seqpos( rsd, rsd_seqpos-1, false );
					jump_pose.insert_residue_by_jump( rsd, rsd_seqpos, rsd_seqpos-1,
																						jump_pose.residue( rsd_seqpos-1 ).atom_index( jump_anchor_atom ), // anchor
																						rsd.atom_index( jump_anchor_atom ) ); // root_atomno
				}

				std::cout << "a " << j << pose.fold_tree() << std::endl;
				pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, rsd.atom_index( jump_anchor_atom ),
																		 partner_rsd.atom_index( jump_anchor_atom ) );
				jump_pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, rsd.atom_index( jump_anchor_atom ),
																					partner_rsd.atom_index( jump_anchor_atom ) );
				std::cout << "b " << j << pose.fold_tree() << std::endl;
			}
			break;
		}
	}

	dump_pdb( pose, "mini_pose.pdb" );
	std::cout << "posefinale: " << pose.fold_tree() << std::endl;
	std::cout << "jump_posefinal: " << jump_pose.fold_tree() << std::endl;
	set_base_partner( pose );
	set_base_partner( jump_pose );

	Pose start_pose, start_jump_pose;
	start_pose = pose;
	start_jump_pose = jump_pose;

	MoveMap mm;
	mm.set_bb( true );
	mm.set_jump( true );


	{
		// setup the options
		MinimizerOptions options( "dfpmin", 0.01, true /*use_nblist*/ );
		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bs, 0.5 );
		scorefxn.set_weight( dna_bp, 0.5 );


		AtomTreeMinimizer minimizer;
		minimizer.run( jump_pose, mm, scorefxn, options );
		dump_pdb( jump_pose, "jump_pose_bs_bp.pdb" );

		jump_pose = start_jump_pose;
	}


	{ // dna_bs minimize
		// setup the options
		MinimizerOptions options( "dfpmin", 0.1, true /*use_nblist*/,
															true /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );
		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bs, 0.5 );


		AtomTreeMinimizer minimizer;
		std::cout << "MINTEST: dna_bs" << std::endl;
		pose = start_pose;
		minimizer.run( pose, mm, scorefxn, options );
		dump_pdb( pose, "pose_bs.pdb" );

		std::cout << "MINTEST: dna_bs jump_pose" << std::endl;
		minimizer.run( jump_pose, mm, scorefxn, options );
		dump_pdb( jump_pose, "jump_pose_bs.pdb" );
	}
	{
		// dna_bs minimize
		MinimizerOptions options( "dfpmin", 0.1, true /*use_nblist*/,
															true /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );
		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bp, 0.5 );


		AtomTreeMinimizer minimizer;
		std::cout << "MINTEST: dna_bp" << std::endl;
		pose = start_pose;
		minimizer.run( pose, mm, scorefxn, options );
		dump_pdb( pose, "pose_bp.pdb" );
	}

}




///////////////////////////////////////////////////////////////////////////////
void
dna_deriv_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace optimization;
	using namespace options;
	using namespace pose;

	using namespace options::OptionKeys;
	using namespace io::pdb;
	using namespace protocols::frags;
	using namespace protocols::moves;
	using namespace scoring::dna;

	using core::Real;

	//
	Pose pdb_pose;
	pose_from_pdb( pdb_pose, "input/1aay.pdb" ); // use 1qn4 for massive dna bend
	Size const nres( pdb_pose.total_residue() );

	set_base_partner( pdb_pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( pdb_pose ) );

	// create a mini-pose with just dna
	Pose pose, jump_pose;

	std::string const jump_anchor_atom( " C2 " );
	for ( Size chain1_begin=1; chain1_begin<= nres; ++chain1_begin ) {
		if ( partner[ chain1_begin ] ) {
			// found first dna residue
			int const chain1( pdb_pose.chain( chain1_begin ) );

			for ( Size j=chain1_begin; ( j<= nres && partner[j] && pdb_pose.chain(j) == chain1 ); ++j ) {
				Residue const &         rsd( pdb_pose.residue(            j ) );
				Residue const & partner_rsd( pdb_pose.residue( partner[ j ] ) );

				Size const rsd_seqpos( j - chain1_begin + 1 );
				if ( rsd_seqpos == 1 ) {
					pose.append_residue_by_bond( rsd );
					jump_pose.append_residue_by_bond( rsd );
				} else {
					pose.append_polymer_residue_after_seqpos( rsd, rsd_seqpos-1, false );
					jump_pose.insert_residue_by_jump( rsd, rsd_seqpos, rsd_seqpos-1,
																						jump_pose.residue( rsd_seqpos-1 ).atom_index( jump_anchor_atom ), // anchor
																						rsd.atom_index( jump_anchor_atom ) ); // root_atomno
				}

				std::cout << "a " << j << pose.fold_tree() << std::endl;
				pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, rsd.atom_index( jump_anchor_atom ),
																		 partner_rsd.atom_index( jump_anchor_atom ) );
				jump_pose.insert_residue_by_jump( partner_rsd, rsd_seqpos + 1, rsd_seqpos, rsd.atom_index( jump_anchor_atom ),
																					partner_rsd.atom_index( jump_anchor_atom ) );
				std::cout << "b " << j << pose.fold_tree() << std::endl;
			}
			break;
		}
	}

	std::cout << "posefinale: " << pose.fold_tree() << std::endl;
	std::cout << "jump_posefinal: " << jump_pose.fold_tree() << std::endl;
	set_base_partner( pose );
	set_base_partner( jump_pose );

	Pose start_pose, start_jump_pose;
	start_pose = pose;
	start_jump_pose = jump_pose;

	MoveMapOP mm ( new MoveMap );
	mm->set_bb( true );
	mm->set_jump( true );

	// setup the options
	scoring::ScoreFunctionOP scorefxn ( new scoring::ScoreFunction );
	moves::MinMover min_mover( mm, scorefxn, "dfpmin", 0.01, true /*use_nblist*/ );

	{
		scorefxn->set_weight( dna_bs, 0.5 );
		scorefxn->set_weight( dna_bp, 0.5 );

		min_mover.apply( jump_pose );
		jump_pose = start_jump_pose;
	}


	{ // dna_bs minimize
		// setup the options
		min_mover.min_type( "linmin" );
		min_mover.tolerance( 0.1 );
		min_mover.deriv_check( true );

		scorefxn->reset();
		scorefxn->set_weight( dna_bs, 0.5 );

		std::cout << "MINTEST: dna_bs" << std::endl;
		pose = start_pose;
		min_mover.apply( pose );

		std::cout << "MINTEST: dna_bs jump_pose" << std::endl;
		min_mover.apply( jump_pose );
	}
	{
		// dna_bs minimize
		scorefxn->reset();
		scorefxn->set_weight( dna_bp, 0.5 );

		std::cout << "MINTEST: dna_bp" << std::endl;
		pose = start_pose;
		min_mover.apply( pose );
	}

}




///////////////////////////////////////////////////////////////////////////////
void
dna_io_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace kinematics;
	using namespace options;
	using namespace pose;

	using namespace options::OptionKeys;
	using namespace io::pdb;
	using namespace protocols::frags;
	using namespace protocols::moves;

	using core::Real;

	// read the filenames
	utility::vector1< std::string > files( start_files() );


	// read the files
	util::prof_reset();
	for ( Size n=1; n<= files.size(); ++n ) {
		Pose pose;
		std::cout << "Reading file: " << files[n] << std::endl;
		pose_from_pdb( pose, files[n] );
		for ( Size i=1; i<= pose.conformation().num_chains(); ++i ) {
			std::cout << "chain_sequence: " << files[n] << ' ' << i << ' ' << pose.chain_sequence(i) << std::endl;
		}
		util::prof_show();
	}


}


///////////////////////////////////////////////////////////////////////////////
void
old_simple_conformation_test()
{


	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

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
simple_centroid_test()
{
	using namespace chemical;

	// read centroid residue set
	ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	pose::Pose pose;
	//io::pdb::pose_from_pdb( pose, *rsd_set, "input/test_in_cen.pdb" );
	io::pdb::pose_from_pdb( pose, *rsd_set, "input/test_in.pdb" );

	scoring::ScoreFunction scorefxn;
	scorefxn.set_weight( scoring::env, 1.0 );
	scorefxn.set_weight( scoring::pair, 1.0 );
	scorefxn.set_weight( scoring::cbeta, 1.0 );
	scorefxn.set_weight( scoring::vdw, 1.0 );

	scorefxn( pose );

	scorefxn.accumulate_residue_total_energies( pose );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		std::cout << "vdw: " << I(4,i) << F(9,2,pose.energies().onebody_energies(i)[ scoring::vdw ] ) << std::endl;
	}
	std::cout << "vdw_total: " << F(9,2,pose.energies().total_energies()[ scoring::vdw ] ) << std::endl;

	// if the input file was test_in.pdb, we can compare the placement of centroids in test_cen.pdb
	// with those in test_in_cen.pdb, which has rosetta++-placed centroids
	//
	// looks good!
	dump_pdb( pose, "test_cen.pdb" );

}

///////////////////////////////////////////////////////////////////////////////
void
simple_frag_test()
{
	using namespace protocols::frags;

	using namespace core::options;
	using namespace core::options::OptionKeys;


	// read vall:
	VallData vall( option[ phil::vall_file ] );

	{
		SingleResidueTorsionFragmentLibrary lib;

		vall.get_frags( 200, "LLL", "HHH", 1.0, 1.0, true, true, true, lib );
	}


	// read centroid residue set
	chemical::ResidueTypeSetCAP rsd_set( chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, *rsd_set, "input/test_in.pdb" );

	std::string const sequence( pose.sequence() );

	TorsionFragmentLibrary lib;
	Size const nres( pose.total_residue() );
	Size const frag_size(3);

	lib.resize( nres - frag_size + 1 );
	for ( Size i=1; i<= nres-frag_size+1; ++i ) {
		std::string const frag_seq( sequence.substr(i-1,3) );
		vall.get_frags( 200, frag_seq, "---", 1.0, 0.0, false, false, true, lib[i] );
	}


	{ // try building an ideal peptide:
		using namespace pose;
		using namespace chemical;
		using namespace conformation;

		Pose pose;
		for ( Size i=1; i<= 20; ++i ) {
			ResidueTypeCAPs const & rsd_list( rsd_set->aa_map( static_cast<AA>(i) ) /*BAD*/ );
			for ( ResidueTypeCAPs::const_iterator iter=rsd_list.begin(), iter_end= rsd_list.end(); iter!= iter_end; ++iter ) {
				ResidueType const & rsd_type( **iter );
				if ( ( rsd_type.is_lower_terminus() == ( i == 1 ) ) &&
						 ( rsd_type.is_upper_terminus() == ( i == 20 ) ) ) {
					ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );
					pose.append_residue_by_bond( *new_rsd, true );
				}
			}
		}
		for ( Size i=1; i<= 20; ++i ) {
			pose.set_omega(i,180.0);
		}
		io::pdb::dump_pdb( pose, "test_ideal.pdb" );
	}
}


///////////////////////////////////////////////////////////////////////////////
void
dna_design_test_old( pose::Pose & pose )
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;

	using namespace core::scoring::dna;

	Size const nres( pose.total_residue() );

	// for bp,bs stats:
	utility::vector1< Size > partner;
	find_basepairs( pose, partner );

	{ // test scoring
		set_base_partner( pose );

		ScoreFunction scorefxn;
		scorefxn.set_weight( dna_bp, 1.0 );
		scorefxn.set_weight( dna_bs, 1.0 );

		scorefxn( pose );
		return;
	}


	{ // test basepair params, basestep params
		utility::vector1< Real > params;
		for ( Size i=1; i<= nres; ++i ) {
			if ( partner[i] ) {
				get_base_pair_params( pose.residue(i), pose.residue( partner[i] ), params );
				if ( i<nres && i+1 != partner[i] && partner[i+1] && partner[i+1] == partner[i]-1 ) {
					assert( !pose.residue(i).is_upper_terminus() );
					get_base_step_params( pose.residue(i), pose.residue( i+1 ), params );
				}
			}
		}
	}

	return;

	//
}

///////////////////////////////////////////////////////////////////////////////
void
dna_coupled_rotamer_design_test()
{

	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace io::pdb;

	using namespace core::scoring::dna;
	using namespace core::options;
	using namespace options::OptionKeys;

	Pose pose;
	pose_from_pdb( pose, "input/1aay.pdb" );
	Size const nres( pose.total_residue() );

	set_base_partner( pose );
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );


	ScoreFunction scorefxn;

	//if ( option[ phil::soft_rep ] ) scorefxn.set_etable( FA_STANDARD_SOFT );

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
	scorefxn.set_weight( hack_elec  , 1.00 ); // hack electrostatics

	show_intrachain_energies( pose, scorefxn );

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
				//task->nonconst_residue_task( ii ).prevent_repacking();
				assert( task->pack_residue(ii) );
			} else {
				task->nonconst_residue_task( ii ).allow_aa( na_ade );
				task->nonconst_residue_task( ii ).allow_aa( na_thy );
				task->nonconst_residue_task( ii ).allow_aa( na_gua );
				task->nonconst_residue_task( ii ).allow_aa( na_cyt );
				assert( task->design_residue(ii) );
			}
		}

		{ // setup residue couplings
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

		Energy rottrial_score;
		if ( true ) {
			pack::rotamer_trials( pose, scorefxn, task );

			rottrial_score = scorefxn( pose );

			std::cout << "Completed DNA-design rotamer_trials_test() with new score: " << rottrial_score << " vs orig: " <<
				score_orig << std::endl;

			{ // test for mismatches
				WatsonCrickResidueMatcher m;
				for ( Size i=1;i<= nres; ++i ){
					if ( partner[i]>i ) {
						std::cout << i << pose.residue(i).aa() << ' ' << pose.residue(partner[i]).aa() << std::endl;
						assert( m( pose.residue(i), pose.residue(partner[i])) ); // fails if mismatch
					}
				}
			}
		} else {
			rottrial_score = scorefxn( pose );
		}

		// now try packing
		task->or_include_current( true ); // tmp hack
		pack::pack_rotamers( pose, scorefxn, task);
		Energy pack_score = scorefxn( pose );

		std::cout << "Completed DNA_design pack_rotamers_test() with new score: " << pack_score << " vs orig: " <<
			rottrial_score << std::endl;

		show_intrachain_energies( pose, scorefxn );

		{ // test for mismatches
			WatsonCrickResidueMatcher m;
			for ( Size i=1;i<= nres; ++i ){
				if ( partner[i]>i ) {
					std::cout << i << pose.residue(i).aa() << ' ' << pose.residue(partner[i]).aa() << std::endl;
					assert( m( pose.residue(i), pose.residue(partner[i])));
				}
			}
		}

		// these are useful but cost a little time to get
		//scorefxn.accumulate_residue_total_energies( pose );
		//pose.energies().show( std::cout );

	}
}

///////////////////////////////////////////////////////////////////////////////
void
dna_design_test()
{
	using namespace io::pdb;

	// read the filenames
	utility::vector1< std::string > files( options::start_files() );

	// read the files
	for ( Size n=1; n<= files.size(); ++n ) {
		pose::Pose pose;
		std::cout << "Reading file: " << files[n] << std::endl;
		pose_from_pdb( pose, files[n] );
		std::cout << "file total_residue: " << pose.total_residue() << std::endl;
		if ( pose.total_residue() ) {
			dna_design_test_old( pose );
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


	pose_from_pdb( pose, /*residue_set, */ "input/1aay.pdb" );

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

		std::cout << "Completed SDT rotamer_trials_test() with new score: " << rottrial_score << " vs orig: " <<
			score_orig << std::endl;

		// now try packing
		task->or_include_current( true ); // tmp hack
		pack::pack_rotamers( pose, scorefxn, task);
		Energy pack_score = scorefxn( pose );

		std::cout << "Completed SDT pack_rotamers_test() with new score: " << pack_score << " vs orig: " <<
			rottrial_score << std::endl;

		// these are useful but cost a little time to get
		//scorefxn.accumulate_residue_total_energies( pose );
		//pose.energies().show( std::cout );

	}


}

///////////////////////////////////////////////////////////////////////////////
void
simple_copy_test()
{

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	pose::Pose new_pose;
	new_pose = pose;

	dump_pdb( new_pose, "test_copy.pdb" );
}


///////////////////////////////////////////////////////////////////////////////
void
simple_hbond_test()
{
	using namespace optimization;
	using id::AtomID;
	using id::DOF_ID;
	using id::PHI;
	using id::THETA;
	using id::D;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	dump_hbond_pdb( pose, "test_out.hbonds.pdb" );

	// setup scorefxn
	scoring::ScoreFunctionOP scorefxn ( new scoring::ScoreFunction );
	scorefxn->set_weight( scoring::hbond_lr_bb, 1.0 );
	scorefxn->set_weight( scoring::hbond_sr_bb, 1.0 );
	scorefxn->set_weight( scoring::hbond_bb_sc, 1.0 );
	scorefxn->set_weight( scoring::hbond_sc, 1.0 );

	(*scorefxn)( pose );

	//pose.energies().show( std::cout );

	// set the moving dofs
	kinematics::MoveMapOP mm1 ( new kinematics::MoveMap );
	kinematics::MoveMapOP mm2 ( new kinematics::MoveMap );
	kinematics::MoveMapOP mm3 ( new kinematics::MoveMap );
	kinematics::MoveMapOP mm4 ( new kinematics::MoveMap );

	// single backbone
	mm1->set_bb( 4, true );

	// all bb and chi
	mm2->set_bb( true );
	mm2->set_chi( true );

	// single dof
	mm3->set( DOF_ID( AtomID(1,4), PHI ), true );

	// everything!
	mm4->set( PHI, true );
	mm4->set( THETA, true );
	mm4->set( D, true );

	moves::MinMover min_mover( mm1, scorefxn, "dfpmin", 0.001, true /*use_nblist*/,
		true /*deriv_check*/ );

	dump_pdb( pose, "output/before.pdb" );

	//minimizer.run( pose, mm1, scorefxn, options );
	//dump_pdb( pose, "output/after1.pdb" );

	min_mover.movemap( mm2 );
	min_mover.apply( pose );
	dump_pdb( pose, "output/after2.pdb" );

	min_mover.movemap( mm3 );
	min_mover.apply( pose );
	dump_pdb( pose, "output/after3.pdb" );

	min_mover.movemap( mm4 );
	min_mover.apply( pose );
	dump_pdb( pose, "output/after4.pdb" );

}


///////////////////////////////////////////////////////////////////////////////
void
atom_tree_torsion_test()
{
	using id::AtomID;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

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
			id::TorsionType const type( r == 1 ? id::BB :
																					id::CHI );
			Size const n( r == 1 ? rsd.mainchain_atoms().size() : rsd.nchi() );

			for ( Size j=1; j<= n; ++j ) {
				id::TorsionID const tor_id( i, type, j );
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
fa_scorefxn_test()
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
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	std::cout << "pose-score: " << scorefxn( pose ) << std::endl;

	scorefxn.accumulate_residue_total_energies( pose );

}

/// STOP REMOVING FUNCTIONS FROM TEST1
///////////////////////////////////////////////////////////////////////////////
void
rotamer_trials_test()
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
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );
	Energy score_orig = scorefxn( pose );

	//for ( int jj = 1; jj <= pose.total_residue(); ++jj )
	//{
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking();

	clock_t starttime = clock();
	pack::rotamer_trials( pose, scorefxn, task );

	clock_t stoptime = clock();
	std::cout << "rotamer_trials took: " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;

	Energy score = scorefxn( pose );
	std::cout << "Completed rotamer_trials_test() with new score: " << score << " vs orig: " << score_orig << std::endl;

	//pose.energies().show( std::cout );

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
pack_rotamers_test()
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

	{ /* test 1 */
	Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );
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

	scorefxn.set_weight( ref, 1.0 );
	scorefxn.set_weight( p_aa_pp, 0.29 );
	score_orig = scorefxn( pose );

	dump_pdb( pose, "test_packrots.pdb" );
	} /* test1 */

	{ /* test2: design */
		Pose pose;
		io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );
		Energy score_orig = scorefxn( pose );
		pack::task::PackerTaskOP designtask( pack::task::TaskFactory::create_packer_task( pose ));
		designtask->initialize_from_command_line().or_include_current( true );
		designtask->set_bump_check( true );
		utility::vector1< bool > residues_to_repack( pose.total_residue(), false );
		for ( Size ii = 1; ii <= 20; ++ii ) residues_to_repack[ ii ] = true;
		designtask->restrict_to_residues( residues_to_repack );

		clock_t starttime = clock();
		pack::pack_rotamers( pose, scorefxn, designtask);
		clock_t stoptime = clock();
		std::cout << "pack_rotamers with design took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
		Energy design_score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() #2 with new score: " << design_score << " vs orig: " << score_orig << std::endl;
		dump_pdb( pose, "test_design.pdb" );
	} /* test2 */

	{ /* test3: resfile */
		Pose pose;
		io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );
		Energy score_orig = scorefxn( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().read_resfile( "input/test_in.resfile" ).or_include_current( true );
		task->set_bump_check( true );

		clock_t starttime = clock();
		pack::pack_rotamers( pose, scorefxn, task);
		clock_t stoptime = clock();
		std::cout << "pack_rotamers after reading resfile took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << std::endl;
		Energy design_score = scorefxn( pose );

		std::cout << "Completed pack_rotamers_test() #3 with new score: " << design_score << " vs orig: " << score_orig << std::endl;
		dump_pdb( pose, "test_design.pdb" );
	} /* test3 */

}

///////////////////////////////////////////////////////////////////////////////
void
small_min_test()
{
	using namespace protocols::moves;
	using namespace scoring;

	using namespace core::options;
	using namespace core::options::OptionKeys;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, start_file() );

	std::string const outfile_prefix( option[ out::file::o ] );

	core::scoring::ScoreFunctionOP scorefxn = new ScoreFunction;

	// aiming for standard packer weights
	scorefxn->set_weight( fa_atr, 0.80 );
	scorefxn->set_weight( fa_rep, 0.44 );
	scorefxn->set_weight( fa_sol, 0.65 );
	scorefxn->set_weight( fa_pair, 0.49 );
	scorefxn->set_weight( fa_dun, 0.56 );
	scorefxn->set_weight( rama, 0.2 );
	scorefxn->set_weight( hbond_lr_bb, 1.17 );
	scorefxn->set_weight( hbond_sr_bb, 1.17 );
	scorefxn->set_weight( hbond_bb_sc, 1.17 );
	scorefxn->set_weight( hbond_sc   , 1.10 );

	// monte carlo object
	MonteCarloOP mc ( new MonteCarlo( pose, *scorefxn, 0.8 /*temperature*/ ) );

	// the movable dof's
	kinematics::MoveMapOP mm ( new kinematics::MoveMap );
	mm->set_bb ( true );
	mm->set_chi( true );

	// options for minimizer
	MinMoverOP min_mover = new MinMover( mm, scorefxn, "dfpmin", 0.001, true /*use_nblist*/ );

	// packer options
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	PackRotamersMoverOP pack_full_repack ( new PackRotamersMover( scorefxn, *task ) );
	/// @bug accumulate_residue_total_energies has to be called before anything can be done
	(*scorefxn)( pose );
	scorefxn->accumulate_residue_total_energies( pose );
	RotamerTrialsMoverOP pack_rottrial ( new EnergyCutRotamerTrialsMover( scorefxn, *task, mc, 0.01 /*energycut*/ ) );
//	pack_rottrial->setup_rottrial_task( pose, mc, 0.01 /*energycut*/ );

	// setup the move objects
	Size nmoves ( 5 );
	SmallMoverOP small_mover( new SmallMover( mm, 0.8/*temp*/, nmoves ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );


	//dump_pdb( pose, "tmp_start.pdb" );
	{ // initial minimization
		TrialMoverOP min_trial ( new TrialMover( min_mover, mc ) );
		min_trial->apply( pose );
	}

	pose::Pose start_pose;
	start_pose = pose;
	int const inner_cycle ( 30 );
	int const outer_cycle ( option[ phil::nloop ] );

	util::prof_reset();

	CycleMoverOP repack_cycle = new CycleMover;

	SequenceMoverOP main_min_seq = new SequenceMover;
	main_min_seq->add_mover( small_mover );
	main_min_seq->add_mover( pack_rottrial );
	main_min_seq->add_mover( min_mover );

	SequenceMoverOP pack_min_seq = new SequenceMover;
	pack_min_seq->add_mover( pack_full_repack );
	pack_min_seq->add_mover( pack_rottrial );
	pack_min_seq->add_mover( min_mover );

	TrialMoverOP main_min_trial = new TrialMover( main_min_seq, mc );
	TrialMoverOP pack_min_trail = new TrialMover( pack_min_seq, mc );
	RepeatMoverOP main_min_cycle = new RepeatMover( main_min_trial, inner_cycle );
	ProfilerMoverOP profiler = new ProfilerMover;

	SequenceMoverOP full_seq = new SequenceMover;
	full_seq->add_mover( main_min_cycle );
	full_seq->add_mover( pack_min_seq );
	full_seq->add_mover( profiler );

	RepeatMoverOP full_min_cycle = new RepeatMover( full_seq, outer_cycle );

	for ( int n=1; n<= option[ out::nstruct ]; ++n ) {
		pose::Pose relax_pose;
		relax_pose = start_pose;

		full_min_cycle->apply( relax_pose );
		util::prof_show();

		if ( outfile_prefix != "none" ) dump_pdb( relax_pose, outfile_prefix+string_of( n )+".pdb" );

	}

}

	//-----------------------------------------------------------------------------------------A
	//_________________________________________________________________________________________A
	//
	//    |
	//    |
	//    |
	//    |
	//
///////////////////////////////////////////////////////////////////////////////
void
rb_test()
{
	using namespace pose;
	using namespace moves;;
	using namespace scoring;
	using namespace kinematics;
	using namespace conformation;
	using namespace chemical;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	int const nres( pose.total_residue() );
	assert( nres > 40 );

	{ // setup a foldtree
		FoldTree f( nres );
		f.new_jump( 8, 40, 15 );
		f.reorder( 8 );
		pose.fold_tree( f );
	}

	dump_pdb( pose, "tmp1.pdb" );

// 	kinematics::Stub stub1( pose.conformation().upstream_jump_stub(1) ),
// 		stub2( pose.conformation().downstream_jump_stub(1) );
// 	exit(0);

	{ // now add a pseudo residue at the end
		ResidueTypeSetCAP residue_set
			( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		ResidueOP rsd( ResidueFactory::create_residue( residue_set->name_map( "VRT" ) ) );
		pose.append_residue_by_jump( *rsd, 1 );
		dump_pdb( pose, "tmp2.pdb" );

		FoldTree f( pose.fold_tree() );
		f.reorder( f.nres() );
		pose.fold_tree( f );

		dump_pdb( pose, "tmp3.pdb" );

		pose.set_jump( 2, Jump() );

		dump_pdb( pose, "tmp4.pdb" );
	}

	// test scoring
	ScoreFunction scorefxn;

	// aiming for standard packer weights
	scorefxn.set_weight( fa_atr, 0.80 );
	scorefxn.set_weight( fa_rep, 0.44 );
	scorefxn.set_weight( fa_sol, 0.65 );
	scorefxn.set_weight( fa_pair, 0.49 );
	scorefxn.set_weight( fa_dun, 0.56 );
	scorefxn.set_weight( rama, 0.2 );
	if ( true ) {
		scorefxn.set_weight( hbond_lr_bb, 1.17 );
		scorefxn.set_weight( hbond_sr_bb, 1.17 );
		scorefxn.set_weight( hbond_bb_sc, 1.17 );
		scorefxn.set_weight( hbond_sc   , 1.10 );
	}

	scorefxn( pose );
//	pose.energies().show( std::cout );

	// now test the simple rigid-body move stuff
	MoveMap mm;
	mm.set_jump( 1, true );

	Pose start_pose;
	start_pose = pose;

	FoldTree f1( pose.fold_tree() ), f2( pose.fold_tree() );
	f1.reorder( 8 );
	f2.reorder( 40 );

	start_pose.fold_tree( f1 );

	// forward rotation
	moves::RigidBodyPerturbMoverOP rb_mover = new RigidBodyPerturbMover(
			pose, 1 /*jump_num*/, 5.0 /*rot*/, 0.0 /*trans*/ );
	moves::PDBDumpMoverOP dumper= new PDBDumpMover( "tmp_fwd_rotation_" );
	moves::SequenceMoverOP sequencer ( new SequenceMover );
	sequencer->add_mover( rb_mover );
	sequencer->add_mover( dumper );

	moves::RepeatMoverOP cycler = new RepeatMover( sequencer, 10 );
	pose = start_pose;
	cycler->apply( pose );

	// forward translation
	pose = start_pose;
	dumper->name( "tmp_fwd_translation_" );
	rb_mover->trans_magnitude( 1.0 );
	rb_mover->rot_magnitude( 0.0 );
	cycler->apply( pose );

	start_pose.fold_tree( f2 );

	// reverse rotation
	pose = start_pose;
	dumper->name( "tmp_rev_rotation_" );
	rb_mover->trans_magnitude( 0.0 );
	rb_mover->rot_magnitude( 5.0 );
	cycler->apply( pose );

	// reverse translation
	pose = start_pose;
	dumper->name( "tmp_rev_tranlation_" );
	rb_mover->trans_magnitude( 1.0 );
	rb_mover->rot_magnitude( 0.0 );
	cycler->apply( pose );
}


///////////////////////////////////////////////////////////////////////////////
void
ccd_test()
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;
	using namespace chemical;
	using namespace optimization;
	using namespace id;

	Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	int const nres( pose.total_residue() );
	assert( nres > 40 );

	int const cutpoint( 18 );
	{ // setup a foldtree
		kinematics::FoldTree f( nres );
		f.new_jump( 8, 26, cutpoint );
		f.reorder( 8 );
		pose.fold_tree( f );

		// introduce new variant types
		if ( true ) {//false ) {
			make_variant_residue( pose, CUTPOINT_LOWER, cutpoint   );
			make_variant_residue( pose, CUTPOINT_UPPER, cutpoint+1 );
		}
	}

	pose.set_phi( 16, pose.phi(16) + 15.0 );
	pose.set_psi( 16, pose.psi(16) - 15.0 );

	dump_pdb( pose, "tmp_before.pdb" );

	Pose start_pose;
	start_pose = pose;

	kinematics::MoveMap mm;

	// setup moving dofs
	for ( int i=15; i<= 22; ++i ) {
		mm.set_bb ( i, true );
		mm.set( TorsionID( i, BB, 3 ), false ); // omega off
	}

	{ // ccd closure

		Real fwd,bwd,tor_delta, rama_delta;
		protocols::loops::fast_ccd_loop_closure( pose, mm, 15, 22, cutpoint, 100, 0.05, true, 100, 100, 100, 100,
																	fwd, bwd, tor_delta, rama_delta );

		dump_pdb( pose, "tmp_after_ccd.pdb" );

	}

	{ // minimization

		// setup the options
		MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/, true /*deriv_check*/ );

		AtomTreeMinimizer minimizer;
		scoring::ScoreFunction scorefxn;
		scorefxn.set_weight( scoring::chainbreak, 1.0 );

		pose = start_pose;
		minimizer.run( pose, mm, scorefxn, options );

		dump_pdb( pose, "tmp_after_dfpmin.pdb" );
	}


}


///////////////////////////////////////////////////////////////////////////////
void
sasa_test()
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in_w_termini.pdb" );

	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;

	Real const total_sasa( scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, 1.4 ) );

	std::cout << "total_sasa: " << total_sasa << std::endl;
	dump_pdb( pose, "dump_test_in_w_termini.pdb" );
}

///////////////////////////////////////////////////////////////////////////////
// should run to completion
//
void
simple_benchmark()
{

	// Simple_min_test(  );
	// Fri Oct 05 13:11:51 EDT 2007 @758 /Internet Time/
	/// Simple_min_test was moved to test/core/optimization/Minimizer.cxxtest.hh

	simple_rotamer_test(  );

	//for (Size ii = 1; ii <= 100; ++ii )
	rotamer_trials_test(  );
	// Wed Oct 10 10:03:07 EDT 2007 @627 /Internet Time/
	// rotamer_trials_test was moved to test/core/pack/RotamerTrials.cxxtest.hh
	// and put back in Oct 20 at 10:09 AM PDS

	pack_rotamers_test(  );

	fa_scorefxn_test(  );

	test_rama(  );

	// simple_conformation_test(  );
	// Thu Nov 08 08:23:43 EST 2007 @599 /Internet Time/
	// simple_conformation_test(  ); was moved to test/core/conformation/Conformation.cxxtest.hh

	atom_tree_torsion_test(  );

	ccd_test(  );

	sasa_test(  );

	simple_dna_test(  );

	dna_coupled_rotamer_design_test();

	dna_deriv_test();

	test_scorefxn_io();
}


///////////////////////////////////////////////////////////////////////////////
void
mm_pack_test()
{
	using namespace conformation;
	using namespace chemical;
	using namespace scoring;
	using namespace pose;

	// scrfxns
	ScoreFunctionOP score_12( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	ScoreFunctionOP score_mod( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	score_mod->set_weight( fa_dun, 0.00 );
	ScoreFunctionOP score_mm_only( new ScoreFunction );	score_mm_only->set_weight( mm_twist,   1.00 );

	// weighting
	Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );

	Energy ener_12  = (*score_12)(pose);
	Energy ener_mod = (*score_mod)(pose);
	Energy ener_mm  = (*score_mm_only)(pose);

	float mm_weight = (ener_12 - ener_mod)/ener_mm;

	ScoreFunctionOP score_mm( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
	score_mm->set_weight( fa_dun, 0.00 );	score_mm->set_weight( mm_twist, mm_weight );

	// pack
	Pose pose_12;
	io::pdb::pose_from_pdb( pose_12, "input/test_in.pdb" );

	Pose pose_mm;
	io::pdb::pose_from_pdb( pose_mm, "input/test_in.pdb" );

// 	Pose pose_mod;
// 	io::pdb::pose_from_pdb( pose_mod, "input/test_in.pdb" );

	pack::task::PackerTaskOP task_12( pack::task::TaskFactory::create_packer_task( pose_12 ));
	task_12->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	task_12->set_bump_check( true );

	pack::task::PackerTaskOP task_mm( pack::task::TaskFactory::create_packer_task( pose_mm ));
	task_mm->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	task_mm->set_bump_check( true );

// 	pack::task::PackerTaskOP task_mod( pack::task::TaskFactory::create_packer_task( pose_mod ));
// 	task_mm->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
//	task_mm->set_bump_check( true );


	Energy orig_score_12 = (*score_12)(pose_12);
	clock_t starttime_12 = clock();	pack::pack_rotamers( pose_12, *score_12, task_12); clock_t stoptime_12 = clock();
	Energy final_score_12 = (*score_12)(pose_12);

	Energy orig_score_mm = (*score_mm)(pose_mm);
	clock_t starttime_mm = clock();	pack::pack_rotamers( pose_mm, *score_mm, task_mm); clock_t stoptime_mm = clock();
	Energy final_score_mm = (*score_mm)(pose_mm);

// 	Energy orig_score_mod = (*score_mod)(pose_mod);
// 	clock_t starttime_mod = clock(); pack::pack_rotamers( pose_mod, *score_mod, task_mod); clock_t stoptime_mod = clock();
// 	Energy final_score_mod = (*score_mod)(pose_mod);

	//output
	std::cout << "mm_pack_test score 12 orig: " << orig_score_12 << " final score: " << final_score_12 << "in "
						<<  ((double) stoptime_12 - starttime_12)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "mm_pack_test score mm orig: " << orig_score_mm << " final score: " << final_score_mm << "in "
						<<  ((double) stoptime_mm - starttime_mm)/CLOCKS_PER_SEC << " seconds" << std::endl;
// 	std::cout << "mm_pack_test score mod orig: " << orig_score_mod << " final score: " << final_score_mod << "in "
// 						<<  ((double) stoptime_mod - starttime_mod)/CLOCKS_PER_SEC << " seconds" << std::endl;
	std::cout << "MM Weight is: " << mm_weight << std::endl;

	dump_pdb(pose_12, "test_out_12.pdb");
	dump_pdb(pose_mm, "test_out_mm.pdb");
//	dump_pdb(pose_mod, "test_out_mod.pdb");

}

///////////////////////////////////////////////////////////////////////////////
void
start_file_test()
{
	utility::vector1< std::string > files( options::start_files() );
	for ( Size i=1; i<= files.size(); ++i ) {
		std::cout << i << ' ' << files[i] << std::endl;
	}
	std::cout << "start_file(): " << options::start_file() << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
// Phil's test function, a little expanded by Rhiju for centroid readin and
// testing multiple pdbs at once.
void
ss_test()
{

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace chemical;

	using namespace scoring;


	ScoreFunction sfxn;
	sfxn.set_weight( hs_pair, 1.0 );
	sfxn.set_weight( ss_pair, 1.0 );
	sfxn.set_weight( rsigma, 1.0 );

	/// Let's get a few PDBs to test
	utility::vector1< std::string > pdbnames( options::start_files() );

	// Centroid! Affects the neighbor list.
	ResidueTypeSetCAP rsd_set( ChemicalManager::get_instance()->residue_type_set( "centroid" ) );

	for ( Size ii = 1; ii < pdbnames.size(); ++ii ) {
		pose::Pose pose;
		std::string const pdbname = pdbnames[ii];

		io::pdb::pose_from_pdb( pose, *rsd_set, pdbname );

		//	protocols::loops::set_secstruct_from_psipred_ss2( pose ); // uses -in::file::psipred_ss2 <ss-filename>

		std::string dsspname = pdbname;
		dsspname.replace( dsspname.rfind(".pdb", dsspname.length() ), 4, ".dssp" );

		protocols::loops::set_secstruct_from_dssp( pose, dsspname );

		std::cout << sfxn( pose ) << std::endl;

		std::cout << pdbname;
		pose.energies().total_energies().show_if_nonzero_weight( std::cout, sfxn.weights() );
		std::cout << std::endl;
	}

}


///////////////////////////////////////////////////////////////////////////////
void
backrub_min_test()
{
	using namespace id;
	using namespace chemical;
	using namespace conformation;
	using namespace kinematics;

	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, "input/test_short.pdb" );

	utility::vector1< AtomID > mainchain;
	Size const seqpos( 3 );

	mainchain.push_back( AtomID( pose.residue( seqpos-1 ).atom_index("CA"), seqpos-1 ) );
	mainchain.push_back( AtomID( pose.residue( seqpos-1 ).atom_index( "C"), seqpos-1 ) );
	mainchain.push_back( AtomID( pose.residue( seqpos   ).atom_index( "N"), seqpos   ) );
	mainchain.push_back( AtomID( pose.residue( seqpos   ).atom_index("CA"), seqpos   ) );
	mainchain.push_back( AtomID( pose.residue( seqpos   ).atom_index( "C"), seqpos   ) );
	mainchain.push_back( AtomID( pose.residue( seqpos+1 ).atom_index( "N"), seqpos+1 ) );
	mainchain.push_back( AtomID( pose.residue( seqpos+1 ).atom_index("CA"), seqpos+1 ) );

	AtomID const downstream_id( AtomID( pose.residue( seqpos+1 ).atom_index("C"), seqpos+1 ) );

	utility::vector1< std::pair< Size, Size > > edges;
	edges.push_back( std::make_pair( 1, 7 ) ); // CA of seqpos-1 --> CA of seqpos+1
	edges.push_back( std::make_pair( 7, 4 ) ); //       seqpos+1 -->       seqpos
	edges.push_back( std::make_pair( 4, 1 ) ); //       seqpos   -->       seqpos-1

	ResidueTypeSet rsd_set( pose.residue(1).residue_type_set() );
	ResidueOP psd( ResidueFactory::create_residue( rsd_set.name_map( "VRT1" ) ) );
	pose.append_residue_by_jump( *psd, 1 );
	pose.append_residue_by_jump( *psd, 1 );
	pose.append_residue_by_jump( *psd, 1 );

	Size const first_new_pseudo_residue( pose.total_residue() - 2 );


	pose.conformation().setup_backrub_segment( mainchain, downstream_id, edges, first_new_pseudo_residue );

	dump_pose_kinemage( "tmp2.kin", pose );


	pose::Pose start_pose;
	start_pose = pose;
	for ( Size i=1; i<= edges.size(); ++i ) {
		pose = start_pose;
		AtomID const id( 1, first_new_pseudo_residue + i - 1 );
		kinematics::Atom const & atom( pose.atom_tree().atom( id ) );
		assert( atom.n_children() > 0 );

		DOF_ID const dof( atom.child(0)->id(), PHI );
		std::cout << "pseudo 1st child: " << i << ' ' <<dof << std::endl;
		Real const orig( pose.dof( dof ) );
		for ( int k=-5;k<=5; ++k ) {
			pose.set_dof( dof, orig+k*0.05 );
			pose.dump_pdb("backrub_move_"+string_of(i)+"_"+string_of(k)+".pdb" );
		}
	}

	pose = start_pose;

	{ // try minimizing those dofs
		using namespace optimization;
		using namespace scoring;

		MoveMap mm;
		for ( Size i=1; i<= edges.size(); ++i ) {
			AtomID const id( 1, first_new_pseudo_residue + i - 1 );
			kinematics::Atom const & atom( pose.atom_tree().atom( id ) );
			assert( atom.n_children() > 0 );

			DOF_ID const dof( atom.child(0)->id(), id::PHI );
			mm.set( dof, true );
		}

		ScoreFunction scorefxn;
		scorefxn.set_weight( fa_atr, 0.80 );
		scorefxn.set_weight( fa_rep, 0.44 );
		scorefxn.set_weight( fa_sol, 0.65 );

		scorefxn.set_weight( hbond_lr_bb, 1.0 );
		scorefxn.set_weight( hbond_sr_bb, 1.0 );
		scorefxn.set_weight( hbond_bb_sc, 1.0 );
		scorefxn.set_weight( hbond_sc, 1.0 );

		pose.dump_pdb( "before.pdb");
		AtomTreeMinimizer().run( pose, mm, scorefxn, MinimizerOptions( "dfpmin", 0.001, true ) );
		pose.dump_pdb( "after.pdb");

	}

// 	Size const first_new_pseudo_residue( pose.total_residue() + 1 );

// 	kinematics::Atom * root( pose.atom_tree().root()->clone(0) );
// 	AtomPointers old_atom_pointer;
// 	root->update_atom_pointer( old_atom_pointer );

// 	for ( Size i=1; i<= mainchain.size(); ++i ) {
// 		assert( mainchain[i] == old_atom_pointer[mainchain[i]]->id() );
// 		std::cout << mainchain[i] << std::endl;
// 	}

// 	kinematics::Atom * new_root
// 		( setup_backrub_atom_tree( mainchain, downstream_id, old_atom_pointer, edges, first_new_pseudo_residue ) );

// 	new_root->show();

// 	AtomTree at( new_root );

// 	dump_atomtree_kinemage( "tmp.kin", at, pose.conformation() );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::init(argc, argv);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////
	// If you want to add a temporary
	// "run this and exit" behavior to
	// main, add it AFTER the following check
	if ( option[ run::benchmark ] ) {
		simple_benchmark(  );
		exit(0);
	}
	/////////////////////////////////////

	backrub_min_test();
	exit(0);

	start_file_test();
	exit(0);

	ss_test();
	exit(0);

	pack_rotamers_test();
	exit( 0 );

	rb_test();
	exit(0);

	test_gb();
	exit(0);

	//dna_io_test();
	//exit(0);

	//list_dihedrals();
	//mm_library_test();
	//mm_score_test();
	mm_pack_test();
	exit(0);

	// rotamer_trials_test(  );
	// Wed Oct 10 10:03:07 EDT 2007 @627 /Internet Time/
	// rotamer_trials_test was moved to test/core/pack/RotamerTrials.cxxtest.hh
	pack_rotamers_test(  );
	exit(0);

	test_scorefxn_io();
	exit(0);

	dna_coupled_rotamer_design_test();
	exit(0);

	dna_deriv_test();
	exit(0);

	dna_design_test();
	exit(0);

	{ // loops graphics
		protocols::viewer::viewer_main( simple_loop_modeling_test_wrapper );
		exit( 0 );
	}

	simple_loop_modeling_test();
	exit(0);

	//{ // loops graphics
		//protocols::viewer::viewer_main( simple_loop_modeling_test_wrapper );
		//exit( 0 );
		//}

	//simple_loop_modeling_test();
	//exit(0);

	//simple_conformation_test();
	//exit(0);

	//simple_frag_test();
	//exit(0);

	//scoring::ScoringManager::get_instance()->get_EnvPairPotential();
	//simple_centroid_test(  );
	//exit(0);

	//simple_dna_test(  );
	//exit(0);

	//set_fullatom_flag_test();
	//exit(0);

	//ccd_test(  );
	//exit(0);

	//small_min_test(  );
	//exit(0);

	//sasa_test(  );

	//rb_test(  );
	//exit(0);

	//exit(0);

	//patch_test( residue_set );
	//exit(0);

	/// Output demo
// 	std::vector<int> A;  A.push_back(1);  A.push_back(2);  A.push_back(3);  A.push_back(5);
// 	utility::vector1<int> B;  B.push_back(10);  B.push_back(20);  B.push_back(30);  B.push_back(45);
// 	std::map<int, std::string> M;  M[1]="one";  M[2]="two";  M[3]="1+2";

// 	T("Demo") << "vector:" << A << " vector1:" << B << " map:" << M << "\n";

// 	T("Error", 10) << "Some error here!!!\n";
// 	T("core.pose") << "Some core pose message\n";
// 	T("core") << "Some core message\n";

// 	Error() << "Some error test...\n";
// 	Warning() << "Some warning test...\n";

	//exit(0);


	//fa_scorefxn_test(  );
	//exit(0);

	//simple_benchmark(  );
	//exit(0);

	//simple_benchmark(  );
	//exit(0);

	//simple_rotamer_test(  );

	//simple_rotamer_test(  );

	//rotamer_trials_test(  );

	//pack_rotamers_test(  );
	//exit(0);

	//fa_scorefxn_test(  );

	//test_rama(  );

	//simple_conformation_test(  );

	//atom_tree_torsion_test(  );

	//exit(0);


	//simple_hbond_test( atom_set,  );
	//exit(0);

	//test_dunbrack_io(  );


	//simple_copy_test(  );
	//exit(0);

	//test_rama( , atom_set );

	// some tests

	//simple_min_test(  );


	//fa_scorefxn_test(  );

}

// {
// 	if ( false ) {
// 		conformation::ConstResidues const & rsd_list
// 			( residue_set.name3_map( "SER" ) );
// 		std::cout << "SER rsd-list: " << rsd_list.size() << std::endl;
// 		std::cout << "ASN rsd-list: " << residue_set.name3_map("ASN").size() <<
// 			std::endl;
// 		std::cout << "aa_asn rsd-list: " <<
// 			residue_set.aa_map( chemical::aa_asn ).size() <<
// 			std::endl;
// 	}


// 	if ( false ) {
// 		// now try reading a residue file
// 		conformation::Residue* rsd_ptr
// 			( conformation::read_topology_file( "ASN.params", atom_set ) );

// 		// try filling in a mini-pose
// 		pose::Pose pose;
// 		pose.append_residue( *rsd_ptr );

// 		pose.fold_tree( kinematics::FoldTree( 1 ) );

// 		dump_pdb( pose, "junk1.pdb" );

// 		pose.set_chi(1,1,0.0);

// 		dump_pdb( pose, "junk2.pdb" );

// 		pose.set_chi(1,1,60.0);

// 		dump_pdb( pose, "junk3.pdb" );

// 	}
// }
// ///////////////////////////////////////////////////////////////////////////////
// // super-simple pdb reader
// //

// void
// read_pdb(
// 	std::string const & filename,
// 	Strings & resids,
// 	Strings & sequence,
// 	Coords & coords
// 	)
// {
// 	using ObjexxFCL::float_of;
// 	typedef numeric::xyzVector< Real > Vector;

// 	resids.clear();
// 	sequence.clear();
// 	coords.clear();


// 	std::ifstream data( filename.c_str() );
// 	if ( !data ) {
// 		std::cout << "unable to open file!! " << filename << std::endl;
// 		return;
// 	}
// 	std::string line;
// 	while ( getline( data, line ) ) {
// 		if ( line.substr(0,6) == "ATOM  " ) {
// 			// parse the info
// 			std::string const atom_name( line.substr(12,4));
// 			std::string const name3( line.substr(17,3));
// 			std::string const resid( line.substr(22,5));

// 			Vector const xyz
// 				( float_of( line.substr(30,8) ),
// 					float_of( line.substr(38,8) ),
// 					float_of( line.substr(46,8) ) );

// 			if ( coords.find( resid ) == coords.end() ) {
// 				resids.push_back( resid );
// 				sequence.push_back( name3 );
// 			}

// 			coords[ resid ] [ atom_name ] = xyz;
// 		}
// 	}
// }


// ///////////////////////////////////////////////////////////////////////////////
// // super-simple
// //
// void
// pose_from_pdb(
// 							pose::Pose & pose,
// 							chemical::ResidueTypeSet const & residue_set,
// 							std::string const & filename
// 							)
// {
// 	//using namespace core;
// 	using namespace conformation;

// 	typedef numeric::xyzVector< Real > Vector;

// 	// reset current data
// 	pose.clear();

// 	Coords coords;
// 	Strings resids, sequence;

// 	read_pdb( filename, resids, sequence, coords );
// 	int const nres_pdb( resids.size() );
// 	for ( int i=1; i<= nres_pdb; ++i ) {
// 		std::string const pdb_name( sequence[i] );
// 		std::string const resid( resids[i] );
// 		std::map< std::string, Vector > const & xyz
// 			( coords.find( resid )->second );

// 		// which residues match this 3letter name?
// 		ResidueTypeCAPs const & rsd_type_list( residue_set.name3_map( pdb_name ) );
// 		if ( rsd_type_list.empty() ) {
// 			std::cout << "Unrecognized aa: " << pdb_name << '\n';
// 			continue;
// 		}

// 		// look for perfect match:
// 		bool matched( false );
// 		for ( Size j=1; j<= rsd_type_list.size(); ++j ) {
// 			ResidueType const & rsd_type( *(rsd_type_list[j]) );

// 			if ( Size( rsd_type.natoms() ) != xyz.size() ) continue;
// 			int mismatches(0);
// 			for ( Size k=1; k<= xyz.size(); ++k ) {
// 				if ( xyz.count( rsd_type.atom_name(k) ) == 0 ) ++mismatches;
// 			}
// 			if ( mismatches ) continue;

// 			matched = true;

// 			// found a perfect match! fill in the coords
// 			ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );

// 			for ( Size k=1; k<= xyz.size(); ++k ) {
// 				new_rsd->atom(k).xyz( xyz.find( new_rsd->atom_name(k) )->second );
// 			}

// 			pose.append_residue( new_rsd );
// 		} // j=1,rsd_type_list.size()
// 		if ( !matched ) {
// 			std::cout << "Unrecognized residue: " << pdb_name << std::endl;
// 		}

// 	} // i=1,nres_pdb

// 	pose.fold_tree( kinematics::FoldTree( pose.total_residue() ) );
// }


// ///////////////////////////////////////////////////////////////////////////////
// void
// patch_test(
// 	chemical::ResidueTypeSet & residue_set
// )
// {
// 	using namespace conformation;
// 	using namespace chemical;
// 	using namespace pose;

// 	residue_set.apply_patches( io::database::full_name( "patches.txt" ) );

// 	Pose pose;
// 	io::pdb::pose_from_pdb( pose, "input/test_in.pdb" );
// 	Size const nres( pose.total_residue() );

// 	if ( false ) {
// 	for ( Size i=2; i<= nres-1; ++i ) {
// 		Residue const & rsd( pose.residue(i) );

// 		ResidueOP new_rsd( rsd.clone() );

// 		for ( Size j=1; j<= rsd.natoms(); ++j ) {

// 			bool is_chi_atom( false );
// 			for ( Size k=1; k<= rsd.nchi(); ++k ) {
// 				if ( j == rsd.chi_atoms(k)[4] ) is_chi_atom = true;
// 			}
// 			if ( is_chi_atom ) continue;

// 			Vector v( rsd.build_atom_ideal( j, pose.conformation() ) ), old_v( rsd.atom(j).xyz() );
// 			Real const dev( v.distance( old_v ) );
// 			if ( dev > 0.5 ) {
// 				std::cout << i << ' ' << rsd.name() << ' ' << j << ' ' << rsd.atom_name(j) << ' ' <<
// 					F(9,3,v.distance(old_v)) << ' ' << old_v << ' ' << v << std::endl;
// 				new_rsd->atom(j).xyz( v );
// 			}
// 		}

// 		//pose.replace_residue( i, *new_rsd, false );
// 	}
// 	}

// 	// build termini
// 	ResidueSelector sel1, sel2;
// 	sel1.add_line( "AA ASP" );
// 	sel1.add_line( "VARIANT_TYPE LOWER_TERMINUS" );
// 	sel1.add_line( "NOT VARIANT_TYPE UPPER_TERMINUS" );

// 	sel2.add_line( "AA LEU" );
// 	sel2.add_line( "PROPERTY UPPER_TERMINUS"); // alternative to variant_type in this case
// 	sel2.add_line( "NOT PROPERTY LOWER_TERMINUS");

// 	ResidueTypeCAPs m1,m2;
// 	residue_set.select_residues( sel1, m1 );
// 	residue_set.select_residues( sel2, m2 );
// 	std::cout << "sel1: " << m1.size() << ' ' << m1[1]->name() << std::endl;
// 	std::cout << "sel2: " << m2.size() << ' ' << m2[1]->name() << std::endl;

// 	ResidueOP nterm_asp( ResidueFactory::create_residue( *(m1[1]), pose.residue(1), pose.conformation() ) );
// 	ResidueOP cterm_leu( ResidueFactory::create_residue( *(m2[1]), pose.residue(nres), pose.conformation() ) );

// 	pose.replace_residue( 1, *nterm_asp, false );
// 	pose.replace_residue( nres, *cterm_leu, false );
// 	dump_pdb( pose, "tmp.pdb" );

// 	id::AtomID_Mask missing;
// 	id::initialize( missing, pose );
// 	for ( Size i=1; i<= nres; ++i ) {
// 		Residue const & rsd( pose.residue(i) );
// 		for ( Size j=1; j<= rsd.natoms(); ++j ) {
// 			if ( rsd.atom_type(j).is_hydrogen() ) {
// 				missing[ id::AtomID(j,i) ] = true;
// 			}
// 		}
// 	}

// 	pose.conformation().fill_missing_atoms( missing );
// 	dump_pdb( pose, "tmp2.pdb" );
// }


// 	// trim back around insertions:
// 	for ( int r=1; r<=2; ++r ) {
// 		for ( int i=1; i< mapping.size1(); ++i ) {
// 			if ( mapping[i] && mapping[i+1] && mapping[i+1] != mapping[i]+1 ) {
// 				std::cout << "insertion: " << r << ' ' << i << ' ' << mapping[i+1] - mapping[i] - 1 << std::endl;
// 				mapping[  i]=0;
// 				mapping[i+1]=0;
// 			}
// 		}
// 		mapping.reverse();
// 	}

// 	std::cout << "trimmed mapping: " << std::endl;
// 	mapping.show();


// 	// setup the initial fold tree:
// 	// identify contiguous blocks of aligned sequence
// 	// make each the core of a peptide segment
// 	//
// 	utility::vector1< std::pair< int, int > > segments;
// 	{
// 		bool in_segment( false );
// 		int seg_begin(0);
// 		for ( Size i=1; i<= mapping.size1(); ++i ) {
// 			if ( mapping[i] ) {
// 				if ( !in_segment ) {
// 					in_segment = true;
// 					seg_begin = i;
// 				}
// 			} else {
// 				if ( in_segment ) {
// 					in_segment = false;
// 					segments.push_back( std::make_pair( seg_begin, i-1 ) );
// 				}
// 			}
// 		}
// 		if ( in_segment ) segments.push_back( std::make_pair( seg_begin, nres1 ) );
// 	} // scope



// 	kinematics::FoldTree f( pose.total_residue() );
// 	f.new_jump( 8, 26, 18 );
// 	f.reorder( 8 );
// 	pose.fold_tree( f );

// 	dump_pdb( pose, "start.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;

// 	pose.conformation().delete_polymer_residue( 1 );

// 	dump_pdb( pose, "del_0.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;

// 	// prepend an alanine residue
// 	ResidueTypeSetCAP residue_set
// 		( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
// 	ResidueOP ala_rsd( ResidueFactory::create_residue( residue_set->name_map( "ALA" ) ) );

// 	pose.conformation().prepend_polymer_residue_before_seqpos( ala_rsd, 1, true );
// 	pose.set_omega(1,180);

// 	dump_pdb( pose, "ala_1.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;

// 	pose.conformation().prepend_polymer_residue_before_seqpos( ala_rsd, 1, true );
// 	pose.set_omega(1,180);

// 	dump_pdb( pose, "ala_2.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;

// 	pose.conformation().delete_polymer_residue( pose.fold_tree().cutpoint(1) );

// 	dump_pdb( pose, "del_1.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;

// 	pose.conformation().delete_polymer_residue( pose.fold_tree().cutpoint(1) );

// 	dump_pdb( pose, "del_2.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;

// 	pose.conformation().append_polymer_residue_after_seqpos( ala_rsd, pose.fold_tree().cutpoint(1), true );

// 	dump_pdb( pose, "ala_3.pdb" );
// 	std::cout << pose.fold_tree() << std::endl;
