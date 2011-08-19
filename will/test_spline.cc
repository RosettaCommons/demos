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


#include <core/init.hh>
#include <core/types.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/interpolation/spline/CompoundInterpolator.hh>

using utility::vector1;
using core::Real;
using namespace numeric::interpolation::spline;



void
test_simple_interp() {
	vector1<Real> x,y;
	x.push_back( 0.0 );	x.push_back( 1.0 );
	y.push_back( 0.0 );	y.push_back( 1.0 );
	SimpleInterpolator interp( x, y, 0.0, 1.0 );
	for( Real x = 0.1; x < 1.0; x += 0.1 ) {
		Real y,dy;
		interp.interpolate(x,y,dy);
		std::cerr << x << " " << y << " " << dy << std::endl;
	}
}

void
test_compound_interp() {
	CompoundInterpolator cinterp;
	{
		vector1<Real> x,y;
		x.push_back( 0.0 );	x.push_back( 1.0 );
		y.push_back( 0.0 );	y.push_back( 1.0 );
		InterpolatorOP interp1( new SimpleInterpolator( x, y, 0.0, 2.0 ) );
		cinterp.add_range( interp1, 0.0, 1.0 );
	}
	{
		vector1<Real> x,y;
		x.push_back( 1.0 );	x.push_back( 2.0 );
		y.push_back( 1.0 );	y.push_back( 1.0 );
		InterpolatorOP interp1( new SimpleInterpolator( x, y, 2.0, -2.0 ) );
		cinterp.add_range( interp1, 1.0, 2.0 );
	}
	{
		vector1<Real> x,y;
		x.push_back( 2.0 );	x.push_back( 3.0 );
		y.push_back( 1.0 );	y.push_back( 0.0 );
		InterpolatorOP interp1( new SimpleInterpolator( x, y, -2.0, 0.0 ) );
		cinterp.add_range( interp1, 2.0, 3.0 );
	}
	for( Real x = 0.0; x <= 3.0; x += 0.01 ) {
		Real y,dy;
		cinterp.interpolate(x,y,dy);
		std::cerr << x << " " << y << " " << dy << std::endl;
	}

}

void
test_spline_generator() {
	SplineGenerator gen( 0, 0, 0, 3, 0, 0 );
	gen.add_known_value( 2, 1, -2 ); // don't have to be in order
	gen.add_known_value( 1, 1, 2 );
	InterpolatorOP interp = gen.get_interpolator();
	for( Real x = 0.0; x <= 3.0; x += 0.01 ) {
		Real y,dy;
		interp->interpolate(x,y,dy);
		std::cerr << x << " " << y << " " << dy << std::endl;
	}

}

void
test_spline_generator2() {

	SplineGenerator gen( 0, 0.0, 0, 1, 0.0, 0 );

	gen.add_known_value( 0.1, 0.1 );
	gen.add_known_value( 0.31, -0.23 );
	gen.add_known_value( 0.53, 0.4 );
	gen.add_known_value( 0.56, -0.1, 6.0 );

	InterpolatorOP interp = gen.get_interpolator();
	for( Real x = 0.0; x <= 1.0; x += 0.001 ) {
		Real y,dy;
		interp->interpolate(x,y,dy);
		std::cerr << x << " " << y << " " << dy << std::endl;
	}

}


int main (int argc, char *argv[])
{

	core::init( argc, argv );

	// test_simple_interp();

	// test_compound_interp();

	test_spline_generator2();

	return 0;

}
