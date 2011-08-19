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

#include <core/pose/Pose.hh>

#include <core/options/option.hh>

#include <core/util/basic.hh>
#include <core/util/prof.hh>

#include <core/io/database/open.hh>

#include <protocols/moves/BackboneMover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/OutputMovers.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/rigid_body_moves.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/viewer/viewers.hh>

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

/////////////////////////////// TESTING
// #ifdef MAC
// #include <GLUT/glut.h>
// #else
// #include "GL/glut.h"
// #endif

// #include <pthread.h>
/////////////////////////////// TESTING



//silly using/typedef

using namespace core;

using utility::vector1;


///////////////////////////////////////////////////////////////////////////////
void
simple_dna_test()
{
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace kinematics;
	using namespace scoring;
	using namespace optimization;
	using namespace io::pdb;

	Pose pose;


	ResidueTypeSet const & residue_set( *(ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) ) );
	pose_from_pdb( pose, residue_set, "input/1aay.pdb" );
	Size const nres( pose.total_residue() );

	std::cout << pose.sequence() << std::endl;

	// for graphics:
	protocols::viewer::add_conformation_viewer( pose.conformation(), "start_pose" );


	{ // rb+sc minimization
		using namespace protocols::moves;
		using namespace core::options;
		using namespace core::options::OptionKeys;

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
		MonteCarloOP mc( new MonteCarlo( pose, *scorefxn, 0.8 /*temperature*/ ) );

		protocols::viewer::add_monte_carlo_viewer( *mc );

		// the movable dof's
		kinematics::MoveMapOP mm ( new kinematics::MoveMap );
		for ( Size i=1; i<= nres; ++i ) {
			if ( pose.residue(i).is_protein() ) {
				mm->set_bb ( i, true );
				mm->set_chi( i, true );
			}
		}

		// options for minimizer
		MinMoverOP min_mover = new MinMover( mm, scorefxn, "dfpmin", 0.001, true /*use_nblist*/ );

		// packer options
		pack::task::PackerTaskOP task
			( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line(); //.restrict_to_repacking().or_include_current( true );

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose.residue(ii).is_protein() ) {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}
		task->or_include_current( true ); // tmp hack

		PackRotamersMoverOP pack_full_repack ( new PackRotamersMover( scorefxn, *task ) );
		/// @bug accumulate_residue_total_energies has to be called before anything can be done
		(*scorefxn)( pose );
		scorefxn->accumulate_residue_total_energies( pose );
		RotamerTrialsMoverOP pack_rottrial ( new EnergyCutRotamerTrialsMover( scorefxn, *task, mc, 0.01 /*energycut*/ ) );
//		pack_rottrial->setup_rottrial_task( pose, mc, 0.01 /*energycut*/ );

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

		// sequence of movers:
		// small_mover, rotamer trials, minimization (AtomTreeMinimizer)
		SequenceMoverOP main_min_seq = new SequenceMover;
		main_min_seq->add_mover( small_mover );
		main_min_seq->add_mover( pack_rottrial );
		main_min_seq->add_mover( min_mover );

		// combine with mc_boltzmann
		TrialMoverOP main_min_trial = new TrialMover( main_min_seq, mc );

		RepeatMoverOP main_min_cycle = new RepeatMover( main_min_trial, inner_cycle );
		ProfilerMoverOP profiler = new ProfilerMover;

		// sequence of movers:
		// full repack, rotamer trials, minimization (AtomTreeMinimizer)
		SequenceMoverOP pack_min_seq = new SequenceMover;
		pack_min_seq->add_mover( pack_full_repack );
		pack_min_seq->add_mover( pack_rottrial );
		pack_min_seq->add_mover( min_mover );

		TrialMoverOP pack_min_trail = new TrialMover( pack_min_seq, mc );

		SequenceMoverOP full_seq = new SequenceMover;
		full_seq->add_mover( main_min_cycle );
		full_seq->add_mover( pack_min_seq );
		full_seq->add_mover( profiler );

		RepeatMoverOP full_min_cycle = new RepeatMover( full_seq, outer_cycle );

		for ( int n=1; n<= option[ out::nstruct ]; ++n ) {
			pose::Pose relax_pose;
			relax_pose = start_pose;

			protocols::viewer::add_conformation_viewer( relax_pose.conformation(), "relax_pose" );

			full_min_cycle->apply( relax_pose );

			util::prof_show();

			//if ( outfile_prefix != "none" ) dump_pdb( relax_pose, outfile_prefix+string_of( n )+".pdb" );

		}
	}




	{ // try rotamer optimizations //////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	simple_dna_test();

	return 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}
