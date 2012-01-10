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





#include <core/types.hh>





#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>


#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/util.hh>



#include <core/init.hh>

#include <core/io/pdb/pose_io.hh>


#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>


//REMOVE LATER!


// C++ headers

//silly using/typedef


#include <core/util/Tracer.hh>
using core::util::T;
using core::util::Error;
using core::util::Warning;

numeric::random::RandomGenerator RG(15431); // <- Magic number, do not change it!!!

using namespace core;
//using namespace protocols;

//using utility::vector1;

//using io::pdb::dump_pdb;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	//using namespace core::options;
	//using namespace core::options::OptionKeys;

	// setup random numbers and options
	core::init(argc, argv);

	// read the pose
	pose::Pose pose;
	io::pdb::pose_from_pdb( pose, options::start_file() ); // gets filename from -s option

	// this is the numbering system relevant for the resfiles (currently... ie not pdb resnums)
	pose.dump_pdb( "start.pdb" );

	// setup the scorefunction to use
	// could make this commandline configurable
	// not sure if this includes all the terms from rosetta++
	// see file in minirosetta_database/scoring/weights/interface.wts
	//
	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( "interface.wts" ) );


	// perform optimization prior to making the mutation?
	{

		{ // here's an example of minimizing the chi angles of position 3
			using namespace optimization;
			kinematics::MoveMap mm;
			mm.set_chi( 3, true );

			AtomTreeMinimizer().run( pose, mm, (*scorefxn), MinimizerOptions( "dfpmin", 0.001, true ) );
		}

		{ // here's an example of packing with a resfile
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
			task->initialize_from_command_line().read_resfile( "input/prepack.resfile" ).or_include_current( true );

			pack::pack_rotamers( pose, (*scorefxn), task);
		}


		{ // here's an example of packing a subset of residues
			pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
			task->initialize_from_command_line().or_include_current( true );

			// only pack at position 3
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				if ( i == 3 ) {
					task->nonconst_residue_task( i ).restrict_to_repacking();
				} else {
					task->nonconst_residue_task( i ).prevent_repacking();
				}
			}
			pack::pack_rotamers( pose, (*scorefxn), task);
		}

	}

	Real const start_score( (*scorefxn)( pose ) );

	// shows the individual terms, INCLUDING THE SCOREFXN WEIGHTS!!!
	std::cout << "score before mutation: " << start_score << ' ' <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;


	{ // now apply the mutation
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().read_resfile( "input/mutate.resfile" ).or_include_current( true );
		pack::pack_rotamers( pose, (*scorefxn), task);
	}

	{ // perhaps now also apply some post-optimization

		// repack at positions within 6A CB-CB of residue 3
		Vector const & nbr_atom( pose.residue(3).nbr_atom_xyz() );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line().or_include_current( true );

		//
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( nbr_atom.distance( pose.residue(i).nbr_atom_xyz() ) < 6.0 ) {
				task->nonconst_residue_task( i ).restrict_to_repacking();
			} else {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		}
		pack::pack_rotamers( pose, (*scorefxn), task);

	}

	Real const final_score( (*scorefxn)( pose ) );

	std::cout << "score after mutation: " << final_score << ' ' <<
		pose.energies().total_energies().weighted_string_of( scorefxn->weights() ) << std::endl;

	pose.dump_pdb( "final.pdb" );

	exit(0);

}
