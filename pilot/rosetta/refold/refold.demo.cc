// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/rosetta/score/refold.demo.cc
/// @brief  Demo of reading, refolding, and writing a protein
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// libRosetta headers
#include <core/conformation/protein/Protein.hh>
#include <core/conformation/protein/io.hh>
#include <core/conformation/TorsionAngle.hh>
#include <core/conformation/builder/ProteinBuilder.hh>
#include <core/io/pdb/read_pdb.hh>
#include <core/options/option.hh>

// Numeric headers
#include <numeric/conversions.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>


/// @brief Demo of reading, refolding, and writing a protein
int
main( int argc, char * argv[] )
{
	using namespace core;
	using namespace core::conformation;
	using namespace core::conformation::builder;
	using namespace core::conformation::protein;
	using namespace core::io::pdb;
	using namespace core::options;
	using namespace std;
	using numeric::conversions::radians;
	typedef  ProteinBuilder::Proteins  Proteins;

	// Check command line
	if ( argc < 2 ) {
		cerr << "Usage: " << argv[0] << " <protein_file>.pdb" << endl;
		exit( EXIT_FAILURE );
	}

	// Get command line options
	options::initialize().load( argc, argv, true );
	options::process();

	// Read the PDB
	Proteins proteins( read_pdb( argv[1] ) );

	// Grab the first protein
	Protein & protein( **proteins.begin() );

	// Tweak the conformation in 2 places just for the fun of it
	protein.phi( int( protein.size() * 0.25 ) ) += radians( +60.0 ); // +60 degrees
	protein.phi( int( protein.size() * 0.75 ) ) -= radians( -60.0 ); // -60 degrees

	// Refold the protein
	protein.refold_N2C(); // N2C direction fullatom refold

	// Write the PDB
	string const output_file( protein.id() + ".tweaked.pdb" );
	std::ofstream pdb_stream( output_file.c_str() );
	pdb_write( pdb_stream, protein );
	cout << endl << protein.id() << " written to " << output_file << endl;

}
