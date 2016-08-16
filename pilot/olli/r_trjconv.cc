// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file r_trjconv.cc
/// @brief tool to handle and convert pdb/silent rosetta structure output
/// @author Oliver Lange
// libRosetta headers


#include <core/types.hh>
#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/relax_protocols.hh>
#include <protocols/abinitio/FoldConstraints.hh>


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

MY_TRACERS("r_trjconv")

using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;

#include <devel/simple_options/option.hh>
#include <core/scoring/rms_util.hh>
using namespace devel::option;

START_OPT(OptMain,"r_trjconv","Conversion between pdb/silent rosetta outputs")
REG(SOPT("f",input,"silent.out","input file"));
REG(SOPT("o",output,"","silent output file"));
REG(BOPT("fa",full_atom,false,"force full-atom output"));
REG(SOPT("rmsd_target",rmsd_struct,"","compute rmsd to this structure"));
REG(SOPT("tag",tag_selected,"","only extract this tag from silent_file"));
/* REG(BOPT("tags",tag_filter,false,"read tags to be processed from stdin")); */
BEGIN_VAR_LIST
std::string input;
std::string output;
bool full_atom;
std::string rmsd_struct;
std::string tag_selected;
END_OPT(OptMain)
#if 0
class ThisApplication {
public:
	static void register_options();
};

using namespace options;
using namespace options::OptionKeys;


#define OPT(akey)																					\
	core::options::option.add_relevant( akey )

#define NEW_OPT(akey,help,adef)								      		  \
	core::options::option.add( akey , help ).def( adef ); 	\
	OPT( akey )

#define OPT_KEY( type, key )                                    \
namespace core { namespace options { namespace OptionKeys { \
			type##OptionKey const key( #key );		                    \
		} } }

OPT_KEY( Boolean, rerun )
OPT_KEY( Boolean, steal )
OPT_KEY( Boolean, start_extended )
OPT_KEY( File, pca )
OPT_KEY( File, rmsd_target )
OPT_KEY( File, sf )
OPT_KEY( Boolean, viol )
OPT_KEY( Integer, viol_level )
OPT_KEY( String, viol_type )
OPT_KEY( StringVector, tag_selector )

void ThisApplication::register_options() {
}

#endif

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


void* my_main( void * )
{
	using namespace core;
	using namespace options;
	using namespace options::OptionKeys;
	using namespace io::silent;
	using namespace pose;
	OptMain opt;
	ProteinSilentFileData sfd;
	sfd.read_file( opt.input );
	ProteinSilentFileData sfd_out;
	pose::PoseOP rmsd_pose( NULL );
	if ( opt.rmsd_struct.size() ) {
		rmsd_pose = new pose::Pose;
		rmsd_pose->clear();
		io::pdb::pose_from_pdb( *rmsd_pose, opt.rmsd_struct );
	}

	for ( ProteinSilentFileData::const_iterator it=sfd.begin_const(), eit=sfd.end_const(); it!=eit; ++it ) {
		Pose pose;
		std::string tag = it->decoy_tag();
		if ( tag == opt.tag_selected || opt.tag_selected.size()==0 ) {
			std::cerr << tag << std::endl;
			//		it->second->print_conformation( std::cout );
			if ( opt.output.size() > 0) {
				if ( rmsd_pose ) { // calculate RMSD
					pose::Pose pose;
					it->fill_pose( pose );
					core::Real CA_rmsd = core::scoring::CA_rmsd( *rmsd_pose, pose );
					// notify the SilentStruct of its RMSD
					it->add_energy( "CA_rmsd",   CA_rmsd   );
				}
				sfd_out.add_structure ( *it );
			} else {
				it->fill_pose( pose );
				if ( opt.full_atom ) {
					switch_to_residue_type_set( pose, chemical::FA_STANDARD );
				}
				pose.dump_pdb( tag + ".pdb");
			}
		}
	}
	if ( opt.output.size() ) {
		sfd_out.write_all( opt.output );
	};

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

