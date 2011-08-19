// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/rosetta/score/pdb.demo.cc
/// @brief  Demo of building a Protein from a PDB file and writing it out to a PDB file
/// @author Ion Yannopoulos (ion)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// libRosetta headers
#include <core/options/option.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Atoms.hh>
#include <core/conformation/amino/AminoAcid.hh>
#include <core/conformation/amino/AminoAcidFactory.hh>
#include <core/conformation/amino/backbone/keys/BackboneUnitKeys.hh>
#include <core/conformation/amino/keys/AtomKeys.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>



int
main( int /*argc*/, char * argv[] )
{
	using namespace core;
	using namespace core::conformation;
	using utility::vector1;
	using namespace core::pose;

	typedef numeric::xyzVector< Real > Vector;

	// read a pdb file
	std::cout << "Reading file: " << argv[1] << std::endl;


	vector1< std::string > resids;
	vector1< std::string > sequence;

	typedef std::map< std::string, std::map< std::string, Vector > > Coords;
	Coords coords;

	{ // scope
		using ObjexxFCL::float_of;

		std::ifstream data( argv[1] );
		std::string line;
		while ( getline( data, line ) ) {
			if ( line.substr(0,6) == "ATOM  " ) {
				// parse the info
				std::string const atom_name( line.substr(12,4));
				std::string const name3( line.substr(17,3));
				std::string const resid( line.substr(22,5));


				Vector const xyz
					( float_of( line.substr(30,8) ),
						float_of( line.substr(38,8) ),
						float_of( line.substr(46,8) ) );

				if ( coords.find( resid ) == coords.end() ) {
					resids.push_back( resid );
					sequence.push_back( name3 );
				}

				coords[ resid ] [ atom_name ] = xyz;

			}

		}

	}

	Pose pose;


	Size const nres( sequence.size() );
	for ( Size i=1; i<= nres; ++i ) {
		std::string const & resid( resids[i] );
		std::string const & name3( sequence[i] );

		Size const natoms( coords[ resid ].size() );

		// create the aminoacid

		amino::AminoAcidOP aa_ptr
			( amino::AminoAcidFactory::create( name3, amino::BackboneUnitKeys::Mid));
		amino::AminoAcid & aa( *aa_ptr );

		std::cout << resid << ' ' << name3 << ' ' << natoms << ' ' <<
			aa.n_atom() << std::endl;

		std::map< std::string, Vector > const & xyz( coords[ resid ] );

		// fill in the coordinates
		for ( std::map< std::string, Vector >::const_iterator it=xyz.begin(),
						it_end=xyz.end(); it != it_end; ++it ) {
			std::string const & atom_name( it->first );
			assert( amino::AtomKeys::has( atom_name ) );
			amino::AtomKey const & atom_key( amino::AtomKeys::key( atom_name ) );

			// fill in the coords
			aa[ atom_key ].position( it->second );
		}


		std::cout << "test: " << aa.backbone_central_atoms()[1]->number() <<
			std::endl;

		// add to Pose
		pose.append_residue( aa );
	}

	// setup fold_tree
	std::cout << "setting simple fold_tree: " << nres << std::endl;

	kinematics::FoldTree f;
	f.simple_tree( nres );
	pose.fold_tree( f ); // will trigger building the atomtree


}
