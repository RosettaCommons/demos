// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/rosetta/score/score.demo.cc
/// @brief  Demo of scoring a Protein with fullatom_energy
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// libRosetta headers
#include <core/conformation/protein/Protein.hh>
#include <core/conformation/builder/ProteinBuilder.hh>
#include <core/io/pdb/read_pdb.hh>
#include <core/options/option.hh>
#include <core/scoring/fullatom_energy.hh>
#include <core/scoring/initialize.hh>
#include <core/scoring/io/score_file.hh>

// C++ headers
#include <cstdlib>
#include <iostream>
#include <string>


/// @brief Demo of scoring a Protein with fullatom_energy
int
main( int argc, char * argv[] )
{
	using namespace core;
	using namespace core::conformation;
	using namespace core::conformation::builder;
	using namespace core::conformation::protein;
	using namespace core::io::pdb;
	using namespace core::options;
	using namespace core::scoring;
	using namespace core::scoring::io;
	using namespace std;
	typedef  ProteinBuilder::Proteins  Proteins;

	// Check command line
	if ( argc < 3 ) {
		cerr << "Usage: " << argv[0] << " <protein_file>.pdb -database <db_path>" << endl;
		exit( EXIT_FAILURE );
	}

	// Get command line options
	options::initialize().load( argc, argv, true );
	options::process();

	// Scoring initialization
	scoring::initialize();

	// Read the PDB
	Proteins proteins( read_pdb( argv[1] ) );

	// Grab the first protein
	Protein & protein( **proteins.begin() );

	// Score the protein
	cout << endl << "Scoring protein " << protein.id() << endl;
	fullatom_energy( protein );

	// Write the score file
	write_score_file( protein, protein.id() + ".fasc", "native", true );

}
