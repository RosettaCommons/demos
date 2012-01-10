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
/// @author James Thompson

#include "funcfit.hh"

#include <core/init.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/types.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>

#include <numeric/random/random.hh>


// C++ headers
#include <iostream>
#include <string>

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] )
{
	core::init( argc, argv );

	// std::string min_type = "dfpmin"; // linmin, dfpmin, dfpmin_armijo, dfpmin_armijo_nonmonotone
	using namespace core::options;
	using namespace core::options::OptionKeys;

	core::Real min_tol   = option[ james::min_tol ];
	std::string min_type = option[ james::min_type ];
	bool use_nblist      = false; // what is this?
	core::optimization::MinimizerOptions min_options( min_type, min_tol, use_nblist );

	Histogram h( 0.5 );
	for ( int i = 1; i <= 100000; ++i ) {
		// Real num = RG.gaussian() * 10;
		Real num = numeric::random::gaussian() * 10;
		if ( num < 0 ) continue;
		h.add_value( num );
	}

	DensityFit f( h );

	core::optimization::Minimizer minimizer( f, min_options );
	core::optimization::Multivec v( 3 );

	v[1] = 1.0;
	v[2] = 1.0;
	v[3] = 0.5;

	std::cout << "starting values:" << std::endl;
	std::cout << F( 8, 3, v[1] )
						<< F( 8, 3, v[2] )
						<< F( 8, 3, v[3] )
						<< std::endl;
	std::cout << "init = " << f( v ) << std::endl;

	minimizer.run( v );

	std::cout << "ending values:" << std::endl;
	std::cout << F( 8, 3, v[1] )
						<< F( 8, 3, v[2] )
						<< F( 8, 3, v[3] )
						<< std::endl;
	std::cout << "term = " << f( v ) << std::endl;

	int width = 12;
	std::cout << A( width, "bin_idx" )
						<< A( width, "value"   )
						<< A( width, "dens" )
						<< A( width, "fit_dens" )
						<< std::endl;

	for ( int i = 1; i <= h.nbins(); i++ ) {
		std::cout << I( 8, i )
							<< F( 8, 3, i * h.bin_size() )
							<< F( 8, 3, h.get_density(i) )
							<< std::endl;
	}

	return 0;
} // int main( int argc, char * argv [] )
