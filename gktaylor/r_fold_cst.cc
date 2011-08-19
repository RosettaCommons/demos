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
#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/relax_protocols.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/abinitio/ConcoordConstraint.hh>


#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/viewer/viewers.hh>
#include <core/util/Tracer.hh>

MY_TRACERS("r_fold_cst")

using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;

void
make_pose_from_sequence_(
	std::string sequence,
	chemical::ResidueTypeSet const& residue_set,
	pose::Pose& pose
) {
	using namespace chemical;
	// clear all of the old data in the pose
	pose.clear();

	// setup the pose by appending the appropriate residues residues
	for ( Size seqpos = 1; seqpos <= sequence.length(); ++seqpos ) {
		char aa = sequence[seqpos-1]; // string indexing is zero-based!
		AA my_aa = aa_from_oneletter_code( aa );
		ResidueTypeCAPs const & rsd_type_list( residue_set.aa_map( my_aa ) );
		Size best_index = 1;
		ResidueType const & rsd_type( *(rsd_type_list[ best_index ]) );
		conformation::ResidueOP new_rsd( conformation::ResidueFactory::create_residue( rsd_type ) );
		if ( seqpos == 1 ) {
			pose.append_residue_by_jump( *new_rsd, 1 );
		} else {
			pose.append_residue_by_bond( *new_rsd, true );
		}
	} // for seqpos
	// pose.conformation().insert_chain_ending( pose.total_residue() - 1 );		// probably not necessary
} // make_pose_match_sequence_

////////////////////////////////////////////////////////////////////////////////////////////////////////////
///@details the function allows a pose to use a different residue_type_set to represent all its residues,
///such as from fullatom residues to centroid residues, or vice versa. During the switch, corresponding atoms
///will be copied. Redundant atoms will be removed (in case from fullatom to centroid) and missing atoms will be
///built by ideal geometry (in the case from centroid to fullatom).
void
switch_to_residue_type_set(
				 pose::Pose & pose,
				 std::string const & type_set_name
)
{
	using namespace core::chemical;
	using namespace core::conformation;

	// retrieve proper residue_type_set
	ResidueTypeSetCAP target_residue_type_set( ChemicalManager::get_instance()->residue_type_set( type_set_name ) );
	// loop each position and find new type that matches from the new type set
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		Residue const & rsd( pose.residue(i) );
		// in future we may have a conformation using mixed type set, so check this by residue
		std::string const & current_type_set_name ( rsd.type().residue_type_set().database_directory() );
		if ( current_type_set_name.find( type_set_name ) != std::string::npos ) {
			std::cerr << "switch_to_residue_type_set: residue " << i << "already in " << type_set_name
								<< " residue_type_set" << '\n';
			continue;
		}
		// get all residue types with same AA
		ResidueTypeCAPs const & rsd_types( target_residue_type_set->aa_map( rsd.aa() ) );
		ResidueOP new_rsd( 0 );
		// now look for a rsdtype with same variants
		for ( Size j=1; j<= rsd_types.size(); ++j ) {
			ResidueType const & new_rsd_type( *rsd_types[j] );
			if ( rsd.type().variants_match( new_rsd_type ) ) {
	new_rsd = ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );
	break;
			}
		}
		if ( ! new_rsd ) {
			std::cerr << "can not find a residue type that matches the residue " << rsd.name()
		<< "at position " << i << '\n';
			utility_exit_with_message( "switch_to_residue_type_set fails\n" );
		}
		// switch to corresponding residue type in the new set.
		pose.replace_residue( i, *new_rsd, false );
	}
}

using namespace core;
void add_constraints (
	pose::Pose &pose,
	std::string fn,
	std::string viol_type = ""
)
{
	using namespace options;
	using namespace options::OptionKeys;
	using namespace scoring::constraints;

	static bool init( false );

	static ConstraintSetOP cst_set( new ConstraintSet() );

	if ( !init ) {
		utility::io::izstream data( fn.c_str() );
		std::string line;
		if ( !data ) {
			std::cerr << "ERROR:: Unable to open constraints file: "
		<< fn << std::endl;
			std::exit( 1 );
		}

		getline(data,line); // header line
		while( getline( data, line ) ) {
			// line format:
			// name1 res1 name2 res2 distance lb ub type
			std::istringstream line_stream( line );
			// std::cout << "line = " << line << std::endl;

			Size res1, res2;
			std::string name1, name2;
			core::Real dist, lb, ub;
			std::string type;
			line_stream
	>> name1 >> res1
	>> name2 >> res2
	>> dist
	>> lb >> ub
	>> type;
			if ( !viol_type.size() || type == viol_type ) {
				trDebug << "add constraint: " << name1 << " " << name2 << " " << res1 << " " << res2 << " " << dist << "\n";
				id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
				id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );
				using protocols::abinitio::CstFunc;
				using protocols::abinitio::CstConstraint;
				using protocols::abinitio::CstType;
				cst_set->add_constraint(  new CstConstraint( atom1, atom2, new CstFunc( lb, ub, 1.0, type ), HBOND ) );
			}
		} // while ( getline(data,line) )
		init = true;
	} // if ( !init )
	trDebug << "add constraint set to pose \n";
	//  cst_set->show( std::cout );
	static ConstraintSetOP cst_setOP = cst_set->clone();
	//  cst_setOP->show(std::cout );
	pose.constraint_set( cst_setOP );
} // add_constraints



#include <devel/simple_options/option.hh>

using namespace devel::option;

START_OPT(OptMain,"r_fold_cst","rosetta with constraints")
REG(COPT("This tool allows to do abinitio folding and relax with constraints. it can also rerun existing outfiles"));
REG(FNOPT("n",native_fn,"native.pdb","native file to determine sequence"));
REG(FNOPT("frag3",frag3,"aaXXXX03","3mer fragments"));
REG(FNOPT("frag9",frag9,"aaXXXX09","9mer fragments"));
REG(FNOPT("c",constraints_fn,"cnc.cst","file containing constraint information"));
REG(FOPT("fc",cst_weight,1.0,"use this weight for the atom_pair_constraints"));
REG(IOPT("nstruct",nstruct,1,"number of structures to generate"));
REG(FOPT("cycles",cycles,1.0,"factor to increase/decrease cycle number"));
REG(BOPT("abinitio",abinitio,true,"run fold_constraints protocol"));
REG(BOPT("rerun",rerun,false,"dump energy and constraint violations for input structures"));
REG(BOPT("relax",relax,false,"relax structures"));
BEGIN_VAR_LIST
std::string native_fn;
std::string frag3;
std::string frag9;
std::string constraints_fn;
Real cst_weight;
Real cycles;
int nstruct;
bool abinitio;
bool relax;
bool rerun;
END_OPT(OptMain)

START_OPT(OptRerun,"rerun","options for rerun with constraints")
REG(COPT("This tool allows to do abinitio folding and relax with constraints. it can also rerun existing outfiles"));
REG(FNOPT("of",silent_fn,"silent.in","go through these structures"));
REG(BOPT("viol",bViol,true,"show violations"));
REG(IOPT("viol_level",viol_level,31,"how much detail for the violation summary"));
REG(BOPT("score",bScore,true,"show score"));
REG(SOPT("tag",tag,"","work only on this tag"));
REG(SOPT("type",viol_type,"","work only on this type of constraints"));
BEGIN_VAR_LIST
std::string silent_fn;
bool bViol;
int viol_level;
bool bScore;
std::string tag;
std::string viol_type;
END_OPT(OptRerun)

void rerun() {
	OptMain opt;
	OptRerun ropt;
	using namespace core;
	using namespace options;
	using namespace options::OptionKeys;
	using namespace io::silent;
	using namespace pose;

	//read silent file for input
	ProteinSilentFileData sfd;
	sfd.read_file( ropt.silent_fn );
	// add native structure to the list
	pose::Pose native_pose;
	native_pose.clear();
	io::pdb::pose_from_pdb( native_pose, opt.native_fn );
	sfd.add_structure( new ProteinSilentStruct( native_pose, "NATIVE", false ) );  /// !!! VERY BAD IF NOT IDEALIZED !!! */
	// run thru all structures
	for ( ProteinSilentFileData::const_iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		Pose pose;
		std::string tag = it->first;
		if ( ropt.tag.size() == 0 || ropt.tag == tag ) {
			it->second->fill_pose( pose,  chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ));
			if ( tag == "NATIVE" ) pose=native_pose; // replace structure with NATIVE so that we don't suffer from non-idealized stuff
			add_constraints( pose, opt.constraints_fn, ropt.viol_type );
			trInfo << tag << " " ;
			pose.dump_pdb("test.pdb");
			if ( ropt.bViol ) {
				pose.constraint_set()->show_violations(  std::cout, pose, ropt.viol_level );
			};
			if ( ropt.bScore ) {
				scoring::ScoreFunctionOP scorefxn
					= scoring::ScoreFunctionFactory::create_score_function( "score3" );
				scorefxn->set_weight( scoring::atom_pair_constraint, opt.cst_weight );
				( *scorefxn )( pose );
				scorefxn->show( std::cout, pose );
				//		it->second->print_conformation( std::cout );
			};
		}
	}
}

void* my_main( void * )
{
	using namespace core;
	using namespace options;
	using namespace options::OptionKeys;
	OptMain opt;
	if (opt.rerun) {
		rerun();
		exit(1);
	};



	pose::Pose native_pose;
	native_pose.clear();
	io::pdb::pose_from_pdb( native_pose, opt.native_fn );

	std::string sequence = native_pose.sequence(); // must match sequence of fragments! maybe read from fasta instead

	ConstantLengthFragSetOP fragset3mer = new ConstantLengthFragSet( 3 );
	ConstantLengthFragSetOP fragset9mer = new ConstantLengthFragSet( 9 );
	fragset3mer->read_fragment_file( opt.frag3 );
	fragset9mer->read_fragment_file( opt.frag9 );

	pose::Pose extended_pose;
	make_pose_from_sequence_( sequence,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )),
		extended_pose
	);

	std::cout << "abinitio start" << std::endl;

	// make extended chain
	for ( Size pos = 1; pos <= extended_pose.total_residue(); pos++ ) {
		extended_pose.set_phi( pos, -45 );
		extended_pose.set_psi( pos, -45 );
		extended_pose.set_omega( pos, 180 );
	}
	add_constraints( extended_pose, opt.constraints_fn );
	// make a MoveMap
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	for ( Size i = 1; i <= extended_pose.total_residue(); ++i ) {
		movemap->set_bb( true );
	}
	FoldConstraints abinitio_protocol( fragset3mer, fragset9mer, movemap );
	abinitio_protocol.init( extended_pose );
	abinitio_protocol.set_cycles( opt.cycles );
	abinitio_protocol.set_constraint_weight( opt.cst_weight );
	// create list of tags already processed

	using protocols::jobdist::BasicJob;
	using protocols::jobdist::BasicJobOP;
	using protocols::jobdist::PlainSilentFileJobDistributor;
	utility::vector1< BasicJobOP > input_jobs;
	int const nstruct = std::max( 1, opt.nstruct );
	BasicJobOP job = new BasicJob("classic_abinitio_relax", nstruct);
	input_jobs.push_back( job );
	PlainSilentFileJobDistributor< BasicJobOP > jobdist( input_jobs );
	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	jobdist.startup();
	while ( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		std::cout << "Starting " << curr_job->output_tag(curr_nstruct) << " ...\n";

		scoring::ScoreFunctionOP scorefxn
			= scoring::ScoreFunctionFactory::create_score_function( "score3" );
		scorefxn->set_weight( scoring::atom_pair_constraint, opt.cst_weight );

		pose::Pose fold_pose ( extended_pose );

		//    protocols::moves::MonteCarlo& mc = abinitio_protocol.mc();
		//    protocols::viewer::add_monte_carlo_silent_viewer( mc, "test_mc_out", false );

		// scoring!
		int start_time = time(NULL);

		// GKT - very important line
		// take logic from here and put into fold_cm_cst
		// setup score function, then take the following line
		// take lines that setup ab initio protocol object
		abinitio_protocol.apply( fold_pose );

		int end_time   = time(NULL);
		std::cout << "TIMEFORABINITIO: " << end_time - start_time << std::endl;

		bool fullatom = false;
		if ( option[ options::OptionKeys::abinitio::relax ].user() ){
			switch_to_residue_type_set( fold_pose, chemical::FA_STANDARD );
			//relax the structure as well
			fullatom = true;

			// change score function to full atom one !
			scoring::ScoreFunctionOP relax_scorefxn = scoring::getScoreFunction();

			ClassicRelax relax_protocol( relax_scorefxn );
			std::cout << "running relax!\n";
			relax_protocol.apply( fold_pose );
		}

		(*scorefxn)(fold_pose);
		jobdist.dump_pose( curr_job, curr_nstruct, fullatom, fold_pose );
		prev_job = curr_job;
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		std::cout << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (pdb_end_time - pdb_start_time) << " seconds.\n";
		abinitio_protocol.clear_checkpoints();
	}
	jobdist.shutdown();
	return NULL;
}

int
main( int argc, char * argv [] )
{
	CommandLineOptions::process_options(argc, argv );
	init( argc, argv );
	protocols::viewer::viewer_main( my_main );
	return 0;
}

