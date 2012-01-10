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
#include <core/fragments/Fragment.hh>
#include <core/fragments/PSSMProfile.hh>
#include <core/fragments/FragmentSet.hh>
#include <core/fragments/FragmentDatabase.hh>

#include <ObjexxFCL/formatted.io.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

int main (int argc, char const* argv[])
{
	using std::string;
	using namespace core::fragments;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::fmt;

	if ( argc != 4 ) {
		std::cout << "usage: pick_frags vall_file profile size" << std::endl;
		std::cout << " - vall_file is the vall database" << std::endl;
		std::cout << " - profile is the .checkpoint formatted profile" << std::endl;
		std::cout << " - size is an integer for fragment length" << std::endl;
		std::exit(1);
	}

	FragmentDatabase fdb;
	string vall_file( argv[1] );
	string profile( argv[2] ) ;
	int size = int_of( argv[3] );

	PSSMProfile myProfile;
	myProfile.read_profile( profile );
	//std::cout << myProfile << std::endl;

	fdb.read_vall( vall_file );
	FragmentSet fs = fdb.make_fragment_set( myProfile, size, 200 );
	fs.print_fragments();

	return 0;
}
