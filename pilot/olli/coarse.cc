// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>

#include <core/util/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
#include <core/scoring/dunbrack/RotamerLibraryScratchSpace.hh>
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

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>

#include <core/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <protocols/relax_protocols.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
//#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>



// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef
#include <ostream>

#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;



using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;

using utility::vector1;

using io::pdb::dump_pdb;

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

	RotamerLibraryScratchSpace scratch;

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		rot_lib.rotamer_energy( pose.residue(i), scratch );
	}
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
void
simple_conformation_test()
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


void set_my_score_weights(scoring::ScoreFunction &scorefxn, std::string tag) {
	using namespace scoring;
	/* NOTE: the force-matching was done with fa_atr=fa_rep=fa_solv=1.0
				 h-bond forces were ignored

				 repeat force-matching for the canonical weights 0.8/0.44/0.65 and h-bonds

				 also no extra parameters for the termini...
				 I edited the patches to get rid of replacement of termini-atoms by some "unknown" atoms

				 either we ignore termini or we also force match termini
	*/

	/*	scorefxn.set_weight( fa_atr, 0.80 );
		scorefxn.set_weight( fa_rep, 0.44 );
		scorefxn.set_weight( fa_sol, 0.65 );
		scorefxn.set_weight( fa_dun, 0.49 );
	*/
	if ( tag == "fine" ) {
		/*		scorefxn.set_weight( fa_atr, 1.0 );
		scorefxn.set_weight( fa_rep, 1.0 );
		scorefxn.set_weight( fa_sol, 1.0 );
		scorefxn.set_weight( fa_dun, 1.0 );*/
		scorefxn.set_weight( fa_atr, 0.8 );
		scorefxn.set_weight( fa_rep, 0.44);
		scorefxn.set_weight( fa_sol, 0.65);
		scorefxn.set_weight( fa_intra_rep, 0.004);
		scorefxn.set_weight( fa_pair, 0);
		scorefxn.set_weight( fa_dun, 0.56);
		scorefxn.set_weight( hbond_lr_bb, 1.17);
		scorefxn.set_weight( hbond_sr_bb, 1.17);
		scorefxn.set_weight( hbond_bb_sc, 1.17);
		scorefxn.set_weight( hbond_sc, 1.1);
		scorefxn.set_weight( p_aa_pp, 0.64);
	}	else if ( tag == "coarse" ) {
		 scorefxn.set_weight( coarse_fa_atr, 1 );
		 scorefxn.set_weight( coarse_fa_rep, 1 );
		 scorefxn.set_weight( coarse_fa_sol, 1 );
		 scorefxn.set_weight( coarse_beadlj, 1 );
		 scorefxn.set_weight( fa_dun, 1 );
	} else {
		std::cerr << __FILE__<< ' ' << __LINE__ << "unknown energy weight-set" << std::endl;
		exit(1);
	};
}

void coarsify (  pose::Pose &coarse_pose, pose::Pose const& fine_pose ) {
	chemical::ResidueTypeSetCAP coarse_res_set = chemical::ChemicalManager::get_instance()->residue_type_set("coarse_two_bead");
	coarse_res_set->coarsify(coarse_pose, fine_pose);
}

///////////////////////////////////////////////////////////////////////////////
void
coarsify_pose_test( pose::Pose fine_pose )
{
	pose::Pose coarse_pose;
	coarsify( coarse_pose, fine_pose );
	coarse_pose.dump_pdb("coarse.pdb");
}

void simple_score_coarse_test( pose::Pose fine_pose ) {
	pose::Pose coarse_pose;
	coarsify( coarse_pose, fine_pose );

	using namespace scoring;
	ScoreFunction scorefxn;
	ScoreFunction coarse_scfxn;

	set_my_score_weights(scorefxn,"fine");
	set_my_score_weights(coarse_scfxn,"coarse");

	Real normal_score, coarse_score;
	normal_score=scorefxn( fine_pose );

	std::cout << " normal pose-score " << normal_score << std::endl;
	std::cout << " fa_dun: " << fine_pose.energies().total_energies()[ fa_dun ];
	std::cout << " coarse pose-score: " << ( coarse_score=coarse_scfxn(coarse_pose) ) << std::endl;

	std::cout << " fa_atr: " << coarse_pose.energies().total_energies()[ coarse_fa_atr ];
	std::cout << " fa_rep: " << coarse_pose.energies().total_energies()[ coarse_fa_rep ];
	std::cout << " fa_sol: " << coarse_pose.energies().total_energies()[ coarse_fa_sol ];
	std::cout << " fa_dun: " << coarse_pose.energies().total_energies()[ fa_dun ];
	std::cout << " coarse_lj: " << coarse_pose.energies().total_energies()[ coarse_beadlj ] << std::endl;
}

void simple_coarse_pack_test(pose::Pose fine_pose) {

	using namespace chemical;
	pose::Pose coarse_pose;
	coarsify(coarse_pose, fine_pose);

	using namespace scoring;
	ScoreFunction scorefxn;
	ScoreFunction coarse_scorefxn;

	set_my_score_weights(scorefxn,"fine");
	set_my_score_weights(coarse_scorefxn,"coarse");

	pack::task::PackerTaskOP fine_task
		( pack::task::TaskFactory::create_packer_task( fine_pose ));
	fine_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	pack::task::PackerTaskOP coarse_task
		( pack::task::TaskFactory::create_packer_task( coarse_pose ));
	coarse_task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );

	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	Energy normal_score = scorefxn( fine_pose );
	pack::pack_rotamers( fine_pose, scorefxn, fine_task );
	Energy pack_score = scorefxn( fine_pose );
	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	Energy coarse_normal_score = coarse_scorefxn( coarse_pose );
	pack::pack_rotamers( coarse_pose, coarse_scorefxn, coarse_task);
	Energy coarse_pack_score = coarse_scorefxn( coarse_pose );

	std::cout<< " after packing score: " << pack_score << " orig: " <<  normal_score << std::endl;
	fine_pose.dump_pdb("after_fa_repack.pdb");
	pose::Pose new_coarse;
	coarsify(new_coarse, fine_pose);
	new_coarse.dump_pdb("coarsed_after_fa_repack.pdb");
	coarse_pose.dump_pdb("coarse_repack.pdb");
	std::cout<< " after coarse packing score: " << coarse_pack_score << " orig: " << coarse_normal_score << std::endl;

}

void
copy_pose_test( pose::Pose fine_pose )
{
	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	pose::Pose new_pose(fine_pose);
	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
	std::cerr << new_pose.total_residue() << std::endl;
	std::cerr << __FILE__<< ' ' << __LINE__ << std::endl;
}

void
coarse_min_test( pose::Pose fine_pose )
{
	using namespace optimization;
	using namespace chemical;
	using namespace pose;

	pose::Pose coarse_pose;
	coarsify(coarse_pose, fine_pose);

	kinematics::MoveMap mm;

	// setup moving dofs
	for ( Size i=2; i<= fine_pose.total_residue()-1; ++i ) {
		mm.set_bb ( i, true );
		mm.set_chi( i, true );
	}

	// setup the options
	MinimizerOptions options( "linmin", 10.0, true /*use_nblist*/,
		true /*deriv_check*/ );

	Pose start_pose(coarse_pose);
	AtomTreeMinimizer minimizer;
	{
		scoring::ScoreFunction scorefxn;
		set_my_score_weights(scorefxn,"coarse");
		Pose pose;
		pose = start_pose;
		std::cout << "start min... score " << scorefxn(pose) << std::endl;
		pose.dump_pdb("before_min.pdb");
		std::cout << "MINTEST: coarse vdw and fa_dun " << std::endl;
		minimizer.run( pose, mm, scorefxn, options );
		pose.dump_pdb("after_min.pdb");
		std::cout << "finished min.... score " << scorefxn(pose) << std::endl;
	}

}

void coarse_relax_test( pose::Pose fine_pose ) {
	using namespace optimization;
	using namespace chemical;
	using namespace pose;

	pose::Pose coarse_pose;
	coarsify(coarse_pose, fine_pose);

	scoring::ScoreFunctionOP coarse_scorefxn=new scoring::ScoreFunction;
	set_my_score_weights( *coarse_scorefxn,"coarse" );

	protocols::ClassicRelax myrelax( coarse_scorefxn);
	myrelax.use_coarse_vdw();

	std::cerr << "RUNNING RELAX" << std::endl;
	myrelax.apply(coarse_pose );

	std::cerr << "DONE" << std::endl;

}

class MyJob {
public:
	virtual void
	run( std::string filename, pose::Pose & pose ) = 0;

	virtual ~MyJob() {};
};

class RescoreStructures : public MyJob {
public:
	RescoreStructures( std::string);

	virtual
	~RescoreStructures() {};

	virtual
	void
	init( pose::Pose  ) { bFirst = false; };

	void
	run( std::string filename, pose::Pose & pose );

	void
	evaluate( std::string filename, pose::Pose &fine, pose::Pose &coarse, std::string mode );

	void
	set_score_weights( scoring::ScoreFunction &scorefxn, std::string tag );
protected:
	virtual void
	do_stuff( pose::Pose &, scoring::ScoreFunctionOP, std::string ) {};

private:
	scoring::ScoreFunctionOP fine_scorefxn;
	scoring::ScoreFunctionOP coarse_scorefxn;
	utility::io::ozstream out;
	bool bFirst;
};

class MinimizeStructures : public RescoreStructures {
public:
	MinimizeStructures( std::string outn ) :
		RescoreStructures (outn),
		min_options_( "linmin", 10.0, true /*use_nblist*/,true /*deriv_check*/ ) { };
protected:
	void
	init( pose::Pose pose );

	void
	do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string tag );

private:
	kinematics::MoveMap mm_;
	optimization::MinimizerOptions min_options_;

};

class RelaxStructures : public RescoreStructures {
public:
	RelaxStructures( std::string out_name ) : RescoreStructures( out_name ) {};
protected:
	void
	do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string tag );
};

class RepackStructures : public RescoreStructures {
public:
	RepackStructures ( std::string out ) : RescoreStructures(out) {};
protected:
	void
	do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string tag );
};

class StripRelaxStructures : public RelaxStructures {
public:
	StripRelaxStructures(std::string outn ) : RelaxStructures(outn) {};
protected:
	void
	do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string tag );
};

void dump_chi( pose::Pose &pose, std::ostream &out ) {
	for (Size seqpos=1; seqpos<=pose.total_residue(); seqpos++) {
		if (pose.residue_type(seqpos).nchi() > 0) {
			// compute rotamer number from chi
			//pose.residue_type(seqpos).get_RotamerLibrary()->dump_rotinfo(pose.residue(seqpos),out);
			//	out << " || " << seqpos << " " << pose.chi(1,seqpos) ;
			out << std::endl;
		}
	}
}

void dump_chi( pose::Pose &pose, std::string outfile ) {
	utility::io::ozstream out(outfile);
	dump_chi(pose,out);
}

void
MinimizeStructures::init( pose::Pose pose ) {
	// setup moving dofs
	std::cerr << "setup minmap ... " << std::endl;
	for ( Size i=2; i<= pose.total_residue()-1; ++i ) {
		mm_.set_bb ( i, true );
		mm_.set_chi( i, true );
	}
	RescoreStructures::init( pose );
}

void
MinimizeStructures::do_stuff(  pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string  )
{
	optimization::AtomTreeMinimizer minimizer;
	minimizer.run( pose, mm_, *scorefxn, min_options_ );
}

void
RepackStructures::do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string )
{
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	pack::pack_rotamers( pose, *scorefxn, task );
}

void
RelaxStructures::do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string tag )
{
	protocols::ClassicRelax myrelax( scorefxn  );
	if ( tag == "coarse" )  {
		myrelax.use_coarse_vdw();
	};
	std::cerr << "RUNNING RELAX" << std::endl;

	protocols::viewer::add_conformation_viewer( pose.conformation(), "trial_pose" );
	moves::MonteCarloOP mc = myrelax.get_mc( pose );
	protocols::viewer::add_monte_carlo_viewer( *mc );

	myrelax.apply( pose );

	std::cerr << "DONE" << std::endl;
}

void
StripRelaxStructures::do_stuff( pose::Pose & pose, scoring::ScoreFunctionOP scorefxn, std::string tag )
{
	pack::task::PackerTaskOP task
		( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking();

	pack::pack_rotamers( pose, *scorefxn, task );
	RelaxStructures::do_stuff(pose,scorefxn,tag);
}

RescoreStructures::RescoreStructures( std::string out_name ) :
	out ( out_name ),
	bFirst(true)
{
	fine_scorefxn=new scoring::ScoreFunction;
	coarse_scorefxn=new scoring::ScoreFunction;
	set_score_weights( *fine_scorefxn,"fine" );
	set_score_weights( *coarse_scorefxn,"coarse" );

	if (!fine_scorefxn) utility_exit_with_message( "scorefxnOP is NULL" );
	out << "description " << " fine score " << " coarse score" << std::endl;
}

void
RescoreStructures::set_score_weights( scoring::ScoreFunction &scorefxn, std::string tag ) {
	using namespace scoring;
#ifdef SMOOTH_ETABLE_
	if ( tag == "fine" ) {
		//score12
		//weights from standard.wts with corrections from score12.wts_patch
		scorefxn.set_weight( fa_atr, 0.8 );
		scorefxn.set_weight( fa_rep, 0.3054);
		scorefxn.set_weight( fa_sol, 0.5559);
	}	else if ( tag == "coarse" ) {
		scorefxn.set_weight( coarse_fa_atr, 0.8 );
		scorefxn.set_weight( coarse_fa_rep, 0.3054 );
		scorefxn.set_weight( coarse_fa_sol, 0.5559 );
		scorefxn.set_weight( coarse_beadlj, 1 );
	} else {
		std::cerr << __FILE__<< ' ' << __LINE__ << "unknown energy weight-set" << std::endl;
		exit(1);
	};
	//		scorefxn.set_weight( fa_intra_rep, 0.004); //ignore in coarse graining !
	scorefxn.set_weight( fa_pair, -0.3502);
	scorefxn.set_weight( fa_dun, 0.4545);
	scorefxn.set_weight( hbond_lr_bb, 1.075);
	scorefxn.set_weight( hbond_sr_bb, 1.075*0.5); //0.5 from score12 patch
	scorefxn.set_weight( hbond_bb_sc, 1.075);
	scorefxn.set_weight( hbond_sc, 1.075);
	scorefxn.set_weight( p_aa_pp, 0.8803*0.5); //0.5 from score12 patch
	scorefxn.set_weight( rama, 0.2 );  //from score12 patch
#else
	if ( tag == "fine" ) {
		//score12
		//weights from standard.wts with corrections from score12.wts_patch
		scorefxn.set_weight( fa_atr, 0.8 );
		scorefxn.set_weight( fa_rep, 0.44);
		scorefxn.set_weight( fa_sol, 0.65);
	}	else if ( tag == "coarse" ) {
		scorefxn.set_weight( coarse_fa_atr, 0.8 );
		scorefxn.set_weight( coarse_fa_rep, 0.44 );
		scorefxn.set_weight( coarse_fa_sol, 0.65 );
		scorefxn.set_weight( coarse_beadlj, 1 );
	} else {
		std::cerr << __FILE__<< ' ' << __LINE__ << "unknown energy weight-set" << std::endl;
		exit(1);
	};
	//		scorefxn.set_weight( fa_intra_rep, 0.004); //ignore in coarse graining !
	scorefxn.set_weight( fa_pair, 0);
	scorefxn.set_weight( fa_dun, 0.56);
	scorefxn.set_weight( hbond_lr_bb, 1.17);
	scorefxn.set_weight( hbond_sr_bb, 1.17*0.5); //0.5 from score12 patch
	scorefxn.set_weight( hbond_bb_sc, 1.17);
	scorefxn.set_weight( hbond_sc, 1.1);
	scorefxn.set_weight( p_aa_pp, 0.64*0.5); //0.5 from score12 patch
	scorefxn.set_weight( rama, 0.2 );  //from score12 patch
#endif
}

void
RescoreStructures::run(std::string filename, pose::Pose & fine_pose ) {
	if (bFirst)
		init(fine_pose);

	pose::Pose coarse_pose;
	coarsify(coarse_pose, fine_pose);

	evaluate( filename, fine_pose, coarse_pose, "pre" );


	do_stuff( fine_pose, fine_scorefxn, "fine" );
	do_stuff( coarse_pose, coarse_scorefxn, "coarse" );
	pose::Pose coarsed_fine_pose;
	coarsify(coarsed_fine_pose, fine_pose);
	dump_pdb (coarsed_fine_pose, filename + ".fine_crsd.pdb");
	dump_chi (fine_pose,filename + ".fine.chi");
	dump_chi (coarsed_fine_pose, filename + ".fine_crsd.chi");
	dump_chi (coarse_pose, filename + ".coarse.chi");
	evaluate(filename, fine_pose, coarse_pose, "post" );
}

void
RescoreStructures::evaluate(
	std::string filename,
	pose::Pose &fine,
	pose::Pose &coarse,
	std::string mode
)
{
	Real fine_score = (*fine_scorefxn)( fine );
	Real coarse_score = (*coarse_scorefxn)( coarse );
	std::cout << "total_energy summary for fine " << std::endl;
	fine.energies().show(std::cout);
	std::cout << "total_energy summary for coarse " << std::endl;
	coarse.energies().show(std::cout);
	out << filename << " " << fine_score << " " << coarse_score << " ";
	if ( mode == "post" ) {
		out << std::endl;
		dump_pdb( fine , filename + ".fine" +".pdb" );
		dump_pdb( coarse , filename + ".coarse" +".pdb" );
	} else {
		dump_pdb( coarse, filename + ".coarse.pre" + ".pdb");;
	}
}

void process_structure_list( std::string file_list, MyJob &job ) {
	using namespace std;
	utility::io::izstream list( file_list );
	if ( !list.good() ) {
		cerr << "can't read file " << file_list << endl;
		exit(1);
	}
	while ( list ) {
		string pdbfile;
		list >> pdbfile;
		if ( strip_whitespace(pdbfile) != "" ) {
			cerr << "read pdb " << pdbfile << "..." << std::endl;
			pose::Pose pose;
			io::pdb::pose_from_pdb( pose, pdbfile );
			int pos = pdbfile.find_last_of("/");
			std::cerr << "str pos of slash = " << pos << std::endl;
			if (pos>0) {
				job.run(pdbfile.substr(pos), pose);
			} else {
				job.run(pdbfile,pose);
			}
		}
	};

}

///////////////////////////////////////////////////////////////////////////////
/*void*
	my_main( void* ) // int argc, char * argv [] ) */
int
main(int argc, char** argv )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

	core::init(argc, argv);

	if ( option[ run::benchmark ] ) {

		std::string const pdbfile ( "input/test_in.pdb" );
		std::cout << "coarse start" << std::endl;
		pose::Pose pose;
		io::pdb::pose_from_pdb( pose, pdbfile );

		//coarsify_pose_test( pose );


		//simple_score_coarse_test( pose );

		//simple_coarse_pack_test( pose );

		//coarse_min_test( pose );

		coarse_relax_test( pose );
	} else {
		std::string jobname = option[ OptionKeys::coarse::cjob ];
		std::string outname = jobname + ".out";
		RescoreStructures *ajob=NULL;
		if ( option[ OptionKeys::coarse::crescore ] )
			ajob = new RescoreStructures(outname);
		else if ( option[ OptionKeys::coarse::crelax ] )
			ajob = new RelaxStructures(outname);
		else if ( option[ OptionKeys::coarse::crepack ] )
			ajob = new RepackStructures(outname);
		else if ( option[ OptionKeys::coarse::cstrip_relax ] )
			ajob = new StripRelaxStructures(outname);
		process_structure_list( option[ OptionKeys::coarse::cjob ], *ajob );
	}
	return 0;
}

/*
	void*
	my_main( void* )
	{
	numeric::random::RandomGenerator::initializeRandomGenerators(
	111111, numeric::random::_RND_TestRun_, "ran3");
	std::cerr << "GLUT my_main " << std::endl;
	relax_test();
	return 0;
	}*/

/*
int
main( int argc, char * argv [] )
{
	// options, random initialization
	std::cerr<< "Initialize rosetta core... " << std::endl;
	core::init( argc, argv );
	std::cerr<< "GLUT call viewer_main " << std::endl;
	protocols::viewer::viewer_main( my_main );
	return 0;
}
*/



/* TODO List:

Catch errors if a certain energytable is used on the wrong residueTypeSet
like score(coarse_pose) instead of coarse_score(coarse_pose)

Write Dunbrack Lib if not in database and read it in if it is there ...

*/

