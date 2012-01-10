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
#include "core/fragments/Fragment.hh"
#include "core/fragments/FragmentSet.hh"
#include "core/fragments/FragmentDatabase.hh"

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

int main (int argc, char const* argv[])
{
	using namespace core::fragments;

	if ( argc != 2 ) {
		std::cerr << "usage: frag_test fragment_file";
	}

	FragmentSet fdb;
	std::string filename( argv[1] );
	fdb.read_fragment_file( filename );
	fdb.print_fragments();

	return 0;
}
