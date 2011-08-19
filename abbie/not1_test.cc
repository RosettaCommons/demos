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
#include <devel/dna/util_public.hh>
//#include <devel/dna/util.hh>
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

#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
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
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

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
#include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>

//silly using/typedef


#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;

//numeric::random::RandomGenerator RG(12323); // <- Magic number, do not change it!!!

using namespace core;
using namespace protocols;

using utility::vector1;
using std::string;
using std::cout;
using std::endl;
using io::pdb::dump_pdb;


util::Tracer tt( "demo.abbie.not1_test", util::t_info );


///////////////////////////////////////////////////////////////////////////////
// some silly helper routines:
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
vector1< Size >
parse_pdb_pos( pose::Pose const & pose )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	vector1< string > pdb_pos_list( option[ dna::specificity::pdb_pos ]() );
	vector1< Size >  pose_pos_list;

	for ( Size i=1; i<= pdb_pos_list.size(); ++i ) {
		std::string const & resid( pdb_pos_list[i] );
		int pos;
		char chain;
		std::size_t fpos( resid.find(':') );
		if ( fpos != string::npos ) {
			pos = int_of( resid.substr(0,fpos) );
			if ( fpos == resid.size()-1 ) {
				chain = ' ';
			} else {
				chain = resid[ fpos+1 ];
			}
		} else {
			pos = int_of( resid );
			chain = ' ';
		}
		if ( chain == '_' ) chain = ' ';
		pose_pos_list.push_back( pose.pdb_info()->pdb2pose( chain, pos ) );
	}
	return pose_pos_list;
}

//////////////////////////////////////////////////////////////////////////
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
show_residue_hbonds(
										pose::Pose const & pose,
										Size const seqpos
)
{
	using namespace scoring::hbonds;

	HBondSet const & hbond_set( static_cast< HBondSet const & >( pose.energies().data().get( util::HBOND_SET )));

	for ( Size i=1; i<= Size(hbond_set.nhbonds()); ++i ) {
		HBond const & hb( hbond_set.hbond(i) );
		if ( hb.don_res() == int(seqpos) || hb.acc_res() == int(seqpos) ) {

			std::cout << "RSD_HBOND " <<
				I(4,hb.don_res()) << ' ' << pose.residue( hb.don_res() ).atom_name( hb.don_hatm()) <<
				I(4,hb.acc_res()) << ' ' << pose.residue( hb.acc_res() ).atom_name( hb.acc_atm ()) <<
				F(9,3,hb.energy()) << F(9,3,hb.weight()) << ' ' << hbond_set.allow_hbond(i) << std::endl;
		}
	}

}

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


void
extract_pdb_id(
  std::string const file,
	std::string & pdb_id
)
{
	// extract pdb_id
	Size folder_ind( file.find_last_of("/\\") );
	folder_ind = (folder_ind==std::string::npos)? 0 : folder_ind+1;

	Size ext_ind( file.find_last_of(".") );
	if (ext_ind==std::string::npos)
		ext_ind = file.length();

	 pdb_id= file.substr( folder_ind, ext_ind-folder_ind ) ;
}



///////////////////////////////////////////////////////////////////////////////
// accumulates into energymap
void
retrieve_residue_pair_energies(
	scoring::Energies const & energies,
	Size pos1,
	Size pos2,
	bool & are_they_neighbors,
	scoring::EnergyMap & emap
)
{

	scoring::EnergyGraph const & energy_graph( energies.energy_graph() );
	using namespace scoring;

	if ( pos1 > pos2 ) {
		Size const tmp = pos2;
		pos2 = pos1;
		pos1 = tmp;
	}
	assert( pos2 > pos1 );
	are_they_neighbors = false;
	for ( graph::Graph::EdgeListConstIter
					iru  = energy_graph.get_node( pos1 )->const_upper_edge_list_begin(),
					irue = energy_graph.get_node( pos1 )->const_upper_edge_list_end();
				iru != irue; ++iru ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
		if ( Size( edge->get_second_node_ind() ) == pos2 ) {
			emap += edge->energy_map();
			are_they_neighbors = true;
			break;
		}
	}

	// now looking for longrange energies, eg GB

	{ // GEN BORN
		LREnergyContainerCOP lrec = energies.long_range_container( methods::gen_born_lr );
		if ( lrec && !lrec->empty() ) {

			for ( ResidueNeighborConstIteratorOP
							rni = lrec->const_upper_neighbor_iterator_begin( pos1 ),
							rniend = lrec->const_upper_neighbor_iterator_end( pos2 );
						(*rni) != (*rniend); ++(*rni) ) {
				if ( rni->upper_neighbor_id() == pos2 ) {
					assert( rni->energy_computed() );
					rni->retrieve_energy( emap ); // pbmod
					break;
				}
			}
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
/// @details  Dump a pdb file with b-factor = total (weighted) repulsive energy per atom
/// and with monitors between clashing atoms above clash_threshold
///

void
dump_clash_pdb(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	Real const clash_threshold, // for showing monitors
	std::string const & filename
)
{
	using namespace scoring;
	using namespace chemical;
	using namespace scoring::methods;
	using namespace scoring::etable;
	using namespace scoring::etable::count_pair;

	scorefxn( pose );

	Real const fa_rep_weight( scorefxn[ fa_rep ] );

	// cached energies object
	Energies & energies( pose.energies() );

	// the neighbor/energy links
	EnergyGraph & energy_graph( energies.energy_graph() );
	EtableEnergy const etable_energy
		( *ScoringManager::get_instance()->etable( scorefxn.energy_method_options().etable_type() ),
			scorefxn.energy_method_options() );

	// setup the global atom numbering that would be used for pdb output
	id::AtomID_Map< int > atom_number;
	setup_atom_number( pose, atom_number );

	std::ofstream out( filename.c_str() );
	out << "load pdb inline\n";
	//show_rasmol_hbonds( hbond_set, atom_number, out );

	id::AtomID_Map< Real > clash_score;

	id::initialize( clash_score, pose, 0.0 );

	for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
		conformation::Residue const & rsd1( pose.residue( i ) );
		for ( graph::Graph::EdgeListIter
				iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			conformation::Residue const & rsd2( pose.residue( j ) );

			CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

			for ( Size ii=1; ii<= rsd1.natoms(); ++ii ) {
				id::AtomID const id1( ii, i );
				for ( Size jj=1; jj<= rsd2.natoms(); ++jj ) {
					id::AtomID const id2( jj, j );

					Real weight;
					if ( cpfxn->count( ii, jj, weight ) ) {
						EnergyMap emap;
						Real dsq;
						etable_energy.atom_pair_energy( rsd1.atom(ii), rsd2.atom(jj), weight, emap, dsq );
						Real const repE( emap[ fa_rep] * fa_rep_weight );
						clash_score[ id1 ] += 0.5 * repE;
						clash_score[ id2 ] += 0.5 * repE;

						if ( repE > clash_threshold ) {
							out << "monitor " << atom_number[id1] << ' ' << atom_number[id2] << '\n';
						}
					}
				}
			}
		}
	}

	out << "exit\n";
	io::pdb::dump_bfactor_pdb( pose, clash_score, out );
	out.close();
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
rescore_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;

	using namespace scoring;
	using namespace conformation;
	using namespace pose;

	utility::vector1< std::string > weights_files;
	if ( option[ weights_tag ]() != "none" ) {
		weights_files.push_back( option[ weights_tag ] );
	} else {
		read_list_file( option[ weights_tag_list ], weights_files );
	}

	vector1< string > files( start_files() );
	for ( Size n=1; n<= files.size(); ++n ) {
		std::string const & file( files[n] );

		Pose pose;
		io::pdb::pose_from_pdb( pose, file );
		if ( pose.empty() ) continue;
		scoring::dna::set_base_partner( pose );

		std::string pdb_id;
		extract_pdb_id(file, pdb_id);

		pose.dump_pdb(pdb_id+".rescore_pdb" );

		for ( Size j=1; j<= weights_files.size(); ++j ) {

			std::string tag(pdb_id);
			std::ofstream total_total_out( (tag+".total_total").c_str() );
			std::ofstream total_out( (tag+".total").c_str() );
			std::ofstream sc_sc_out( (tag+".sc_sc").c_str() );

			ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( weights_files[ j ] ) );
			EnergyMap const & weights( scorefxn->weights() );

			// Write files' header
			total_out << '#' << std::setw(9) << " ";
			sc_sc_out << '#' << std::setw(19) << " ";
			for ( Size j=1; j<= n_score_types; ++j ) {
				ScoreType const t =  ScoreType(j);
				if ( weights[ t ] != 0.0 ) {
					total_out << std::setw(14) << t;
					sc_sc_out << std::setw(14) << t;
				}
			}

			vector1< ScoreType > types;
			total_out << std::endl <<  "#" << std::setw(9) << " ";
			sc_sc_out << std::endl <<  "#" << std::setw(19) << " ";
			for ( Size j=1; j<= n_score_types; ++j ) {
				ScoreType const t =  ScoreType(j);
				if ( weights[ t ] != 0.0 ) {
					total_out << F(14, 3, weights[ t ]);
					sc_sc_out << F(14, 3, weights[ t ]);
					types.push_back(t);
				}
			}

			total_out << std::endl;
			sc_sc_out << std::endl ;

			(*scorefxn)( pose );
			Energies const & energies( pose.energies() );
			scorefxn->accumulate_residue_total_energies( pose );

			scorefxn->show( total_total_out, pose );


			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				conformation::Residue const & res( pose.residue( i ) );

				EnergyMap const & emap( energies.residue_total_energies(i) );
				// EnergyMap const & emap( energies.onebody_energies( i ) );

				total_out << std::setw(4) << i << std::setw(2) << res.chain() << std::setw(4) << res.name3();

				vector1< ScoreType >::const_iterator iter;
				for ( iter=types.begin(); iter!=types.end(); ++iter ) {
					ScoreType const t = (*iter);
					total_out << F(14,3,emap[ t ]);
				}
				total_out << std::endl;

				for ( Size i2=i+1; i2<=pose.total_residue(); ++i2 ) {

					conformation::Residue const & res2( pose.residue( i2 ) );

					bool are_they_neighbors;
					EnergyMap pair_emap;

					retrieve_residue_pair_energies( energies, i, i2, are_they_neighbors, pair_emap);

					if ( are_they_neighbors ) {
						sc_sc_out << std::setw(4) << i  << std::setw(2) <<  res.chain() << std::setw(4) <<  res.name3()
											<< std::setw(4) << i2 << std::setw(2) << res2.chain() << std::setw(4) << res2.name3();

						for ( iter=types.begin(); iter!=types.end(); ++iter ) {
							ScoreType const t =  (*iter);
							sc_sc_out << F(14,3,pair_emap[ t ]);
						}
						sc_sc_out << std::endl;
					}
				} // pairwise energies
			}   // sidechain - bb energies
		}
	}
}



///////////////////////////////////////////////////////////////////////////////
void
not1_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;
 	using namespace scoring;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace optimization;
 	using namespace scoring::dna;
 	using namespace pose;



	// read structure

	Pose pose;
	io::pdb::pose_from_pdb( pose, start_file() );


	ScoreFunctionOP scorefxn;
	if ( option[ score_function ].present() ) {
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

	(*scorefxn)(pose);

	{ // minimize sidechains
		scorefxn->show( std::cout, pose );

		Pose min_pose( pose );
		MoveMap mm;
		for ( Size i=1; i<= pose.total_residue(); ++i ) mm.set_chi( i, pose.residue(i).is_protein() );
		AtomTreeMinimizer().run( min_pose, mm, *scorefxn, MinimizerOptions("dfpmin",0.00001,true) );
		scorefxn->show( std::cout, min_pose );
		min_pose.dump_pdb( "minimized.pdb" );

	}


	dump_hbond_pdb( pose, "hbond1.pdb" );
	dump_clash_pdb( pose, *scorefxn, 2.0, "clash1.pdb" );

	{ // try making the modification
		//
		Size const mcy_pos( option[ motif_begin ] );
		assert( pose.residue( mcy_pos ).name() == "CYT" );

		ResidueTypeSet const & residue_set( pose.residue(1).residue_type_set() );
		replace_pose_residue_copying_existing_coordinates( pose, mcy_pos, residue_set.name_map("5MC") );
		pose.dump_pdb("mcy"+string_of( mcy_pos )+".pdb");
		dump_hbond_pdb( pose, "hbond2.pdb" );
		dump_clash_pdb( pose, *scorefxn, 2.0, "clash2.pdb" );

		exit(0);
	}


	(*scorefxn)(pose); // HACK

	{ // try making the thing we want
		pose.set_chi( 2, 178, 180 + pose.chi( 2, 178 ) );
		pose.dump_pdb( "premin.pdb" );
		MoveMap mm;
		mm.set_chi( 178, true );
		AtomTreeMinimizer().run( pose, mm, *scorefxn, MinimizerOptions("dfpmin",0.001,true) );
		pose.dump_pdb( "postmin.pdb" );
	}


	{ // here's an example of packing with a resfile
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().read_resfile( option[ packing::resfile ]() ).or_include_current( true );
		task->set_bump_check( false );
		pack::pack_rotamers( pose, (*scorefxn), task);
		std::cout << "packed score: " << (*scorefxn)( pose ) << std::endl;
		std::cout << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;
		pose.dump_pdb( "packed.pdb" );
	}

	{ // here's an example of packing with a resfile
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().read_resfile( option[ packing::resfile ]() ).or_include_current( true );
		task->set_bump_check( false );
		pack::rotamer_trials( pose, (*scorefxn), task);
		std::cout << "packed_rottrialed score: " << (*scorefxn)( pose ) << std::endl;
		std::cout << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;
		dump_hbond_pdb( pose, "packed_and_rottrialed.pdb" );
		show_residue_hbonds( pose, 178 );

		{
			Residue const & rsd1( pose.residue( 178 ) );
			Residue const & rsd2( pose.residue( 217 ) );

			hbonds::HBondSet hbond_set;
			hbond_set.use_hb_env_dep(true);
			hbond_set.smooth_hb_env_dep(true);

			// prior to finding these bonds, assume rsd1's backbone is
			// donating an hbond rsd2's backbone is accepting an hbond
			// therefore do not identify any hbonds between a backbone and
			// a sidechain that would make a 3-center bond with one of
			// pre-existing hbonds.  Note these restrictions are all
			// protein centric.

			// exclude all backbone-backbone hbonds

			// rsd1 as donor rsd2 as acceptor.  Exclude hbonds where the
			// accepting atom on rsd2 is side chain and donating atom of
			// rsd2 is backbone.
			hbonds::identify_hbonds_1way( rsd1, rsd2, 10, 10, false /*evaluate_derivative*/,
																		true, false, true, false, hbond_set);
			// rsd2 as donor rsd1 as acceptor.  Exclude hbonds where the
			// accepting atom on rsd2 is backbone and the donating atom on
			// rsd1 is sidechain.
			hbonds::identify_hbonds_1way( rsd2, rsd1, 10, 10, false /*evaluate_derivative*/,
																		true, true, false, false, hbond_set);
			TwoBodyEnergyMap hbond_emap;
			get_hbond_energies_new( hbond_set, hbond_emap);

			std::cout << "HBE: "
								<< hbond_emap[ hbond_sc ]	<< ' '
								<< hbond_emap[ hbond_sr_bb ] << ' '
								<< hbond_emap[ hbond_lr_bb ] << ' '
								<< hbond_emap[ hbond_sr_bb_sc ] + hbond_emap[ hbond_lr_bb_sc ] << std::endl;
		}
	}

	{ // replace
		using namespace conformation;
		using namespace chemical;
		ResidueTypeSet const & residue_set( pose.residue(1).residue_type_set() );
		ResidueType const & HIS_E( residue_set.name_map( "HIS" ) );
		Residue const old_rsd( pose.residue( 178 ) );
		assert( old_rsd.name() == "HIS_D" );
		std::cout << old_rsd.name() << '\n';
		ResidueOP new_rsd( ResidueFactory::create_residue( HIS_E, old_rsd, pose.conformation() ) );
		new_rsd->set_chi( 1, old_rsd.chi(1) );
		new_rsd->set_chi( 2, old_rsd.chi(2) );
		pose.replace_residue( 178, *new_rsd, false );
		std::cout << "his_e: " << (*scorefxn)( pose ) << std::endl;
		std::cout << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;
		show_residue_hbonds( pose, 178 );
		dump_hbond_pdb( pose, "his_e.pdb" );
	}

}

///////////////////////////////////////////////////////////////////////////////
void
design_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::options::OptionKeys::dna::specificity;
 	using namespace scoring;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace optimization;
 	using namespace scoring::dna;
 	using namespace pose;

	// read structure

	Pose pose;
	io::pdb::pose_from_pdb( pose, start_file() );


	ScoreFunctionOP scorefxn;
	if ( option[ score_function ].present() ) {
		scorefxn = new ScoreFunction();
		scorefxn->initialize_from_file( option[ score_function ] );
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( "dna.wts" );
	}

	util::prof_reset();

	{ // score the starting structure
		(*scorefxn)(pose);
		std::cout << "start score: " << (*scorefxn)( pose ) << std::endl;
		std::cout << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;
		util::prof_show();
	}


	{ // here's an example of packing with a resfile
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().read_resfile( option[ packing::resfile ]() ).or_include_current( true );
		task->set_bump_check( true );
		pack::pack_rotamers( pose, (*scorefxn), task);

		std::cout << "packed score: " << (*scorefxn)( pose ) << std::endl;
		std::cout << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;

	}

	// write the final pdb
	pose.dump_pdb( option[ out::file::o ] );

	util::prof_show();
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
methylate_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	pose::Pose pose;

	// read the pdb
	io::pdb::pose_from_pdb( pose, start_file() );

	// get the list of positions to be methylated
	vector1< Size > pos_list( parse_pdb_pos( pose ) );

	chemical::ResidueTypeSet const & residue_set( pose.residue(1).residue_type_set() );
	for ( Size j=1; j<= pos_list.size(); ++j ) {
		replace_pose_residue_copying_existing_coordinates( pose, pos_list[j], residue_set.name_map("5MC") );
	}

	pose.dump_pdb( option[ out::file::o ] );


}


///////////////////////////////////////////////////////////////////////////////
void
interface_design_test()
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	if ( !option[ out::file::o ].present() ||
			 !option[ dna::specificity::score_function ].present() ||
			 !option[ dna::specificity::pdb_pos ].present() ||
			 !option[ in::file::s ].present() ) {
		utility_exit_with_message( "Use -s for starting strx, -o for output structure, -pdb_pos for dna positions to be designed around, -score_function for the score weights file");
	}

	// setup scorefxn
	scoring::ScoreFunctionOP scorefxn = new scoring::ScoreFunction();
	scorefxn->initialize_from_file( option[ dna::specificity::score_function ] );


	pose::Pose pose;

	// read the pdb
	io::pdb::pose_from_pdb( pose, start_file() );

	std::cout << "start score: " << (*scorefxn)( pose ) << std::endl;

	// get the list of positions to be methylated
	vector1< Size > dna_pos_list( parse_pdb_pos( pose ) );

	// set up a task to tell the packer what to design
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );
	task->set_bump_check( true );

	Size const nres( pose.total_residue() );
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.is_protein() ) {
			bool close( false ), contact( false );
			devel::dna::check_residue_proximity_to_dna( i, dna_pos_list, pose, close, contact );
			if ( !contact ) {
				if ( close ) {
					task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					task->nonconst_residue_task( i ).prevent_repacking();
				}
			}
		} else {
			// currently not repacking DNA
			task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	{ // debugging
		id::AtomID_Map< Real > bfactor;
		id::initialize( bfactor, pose );

		for ( Size i=1; i<= nres; ++i ) {
			bool const packing( task->pack_residue(i) );
			bool const designing( task->design_residue(i) );
			Real const bf( designing ? 75 : ( packing ? 50 : 25 ) );
			for ( Size j=1; j<= pose.residue(i).natoms(); ++j ) {
				bfactor[ id::AtomID( j, i ) ] = bf;
			}
		}
		std::ofstream out( "pack_or_design.pdb" );
		io::pdb::dump_bfactor_pdb( pose, bfactor, out );
		out.close();
	}

	vector1< std::pair< Real, std::string > > results;
	pack::pack_rotamers_loop( pose, (*scorefxn), task, 50, results );

	std::cout << "redesigned score: " << (*scorefxn)( pose ) << std::endl;
	std::cout << pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;

	string const outfilename( option[ out::file::o ] );
	std::ofstream out( outfilename.c_str() );
	pose.dump_pdb( out );
	scorefxn->show( out, pose );
	out.close();

}


///////////////////////////////////////////////////////////////////////////////
//
// loop over DNA positions specified with -pdb_pos
//
// at each position, try out all 4 base pairs
// for each, rescore, then repack nbring protein residues and rescore
// dump out final pdb and energies
//


void
dna_scan_test()
{



}


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	std::string const mode( option[ dna::specificity::mode ].value() );
	if ( mode == "read_write" ) {
		read_write();
		exit(0);
	} else if ( mode == "not1" ) {
		not1_test();
		exit(0);
	} else if ( mode == "design" ) {
		design_test();
		exit(0);
	} else if ( mode == "rescore" ) {
		rescore_test();
		exit(0);
	} else if ( mode == "dna_scan" ) {
		dna_scan_test();
		exit(0);
	} else if ( mode == "methylate" ) {
		methylate_test();
		exit(0);
	} else if ( mode == "interface_design" ) {
		interface_design_test();
		exit(0);
	}


	exit(0); // add new mode strings

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	// initialize option and random number system
	core::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}

