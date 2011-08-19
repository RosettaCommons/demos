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


#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#include "boinc_zip.h"
#endif // BOINC

#include <core/types.hh>
#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <utility/exit.hh>


std::string DEFAULT_DB = "minirosetta_database_rev19451.zip";

int
main( int argc, char * argv [] )
{
	using namespace core::options;
	using namespace core::options::OptionKeys;

#ifdef BOINC
	utility::boinc::Boinc boinc_wu = utility::boinc::Boinc::instance();
	boinc_wu.initialize( argc, argv );
#endif

	// options, random initialization
	core::init( argc, argv );

#ifdef BOINC
	// unzip a database archive?
	std::string zippedfile = DEFAULT_DB;
	if (option[ in::file::zip ].user())
			zippedfile = option[ in::file::zip ]();
	std::string resolvedfile = zippedfile;
	utility::boinc::resolve_filename( resolvedfile );
	if (!utility::file::file_exists( resolvedfile )) {
		utility_exit_with_message("in::file::zip "+ zippedfile+" does not exist!");
	} else {
		boinc_zip(UNZIP_IT, resolvedfile, "minirosetta_database");
	}
#endif

	if ( option[ run::protocol ]() == "abrelax" )
		return classic_abinitio_relax_main();
	else if ( option[ run::protocol ]() == "ligand_dock" )
		// return ligand_dock_main(argc, argv);
		return ligand_dock_main();
	else {
		utility_exit_with_message("Invalid protocol requested: "+option[ run::protocol ]());
		return 0; // makes compiler happy
	}

}

#ifdef BOINC
#ifdef _WIN32

/*******************************************************
 * Windows: Unix applications begin with main() while Windows applications
 * begin with WinMain, so this just makes WinMain() process the command line
 * and then invoke main()
 */

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst,
									 LPSTR Args, int WinMode)
{
		LPSTR command_line;
		char* argv[1024];
		int argc;

		command_line = GetCommandLine();
		argc = parse_command_line( command_line, argv );
		return main(argc, argv);
}
#endif
#endif
