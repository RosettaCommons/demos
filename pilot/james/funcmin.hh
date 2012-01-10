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

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/types.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/MinimizerOptions.hh>


// C++ headers
#include <iostream>
#include <string>

///////////////////////////////////////////////////////////////////////////////

// f(a,b,c)  = sqrt( (x-a)^2 + (y-b)^2 + (z-c)^2 )
// derivative with respect to a:
// f'(a,b,c) da = (x-a) / f(a,b,c)


class StupidFunc : public core::optimization::Multifunc
{

public:
	StupidFunc( core::Real x, core::Real y, core::Real z ):
	x_( x ),
	y_( y ),
	z_( z )
	{}

	virtual
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const
	{
		core::Real score = std::sqrt(
			( x_ - vars[1] ) * ( x_ - vars[1] ) +
			( y_ - vars[2] ) * ( y_ - vars[2] ) +
			( z_ - vars[3] ) * ( z_ - vars[3] )
		);

		return score;
	}

	virtual
	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const
	{
		using namespace std;
		core::Real func = std::sqrt(
			( x_ - vars[1] ) * ( x_ - vars[1] ) +
			( y_ - vars[2] ) * ( y_ - vars[2] ) +
			( z_ - vars[3] ) * ( z_ - vars[3] )
		);


		dE_dvars[1] = ( vars[1] - x_ ) / func;
		dE_dvars[2] = ( vars[2] - y_ ) / func;
		dE_dvars[3] = ( vars[3] - z_ ) / func;

		// std::cout << "derivatives:" << std::endl
		// 					<< "dE_dvars[1] = " << dE_dvars[1]
		// 					<< "dE_dvars[2] = " << dE_dvars[2]
		// 					<< "dE_dvars[3] = " << dE_dvars[3]
		// 					<< std::endl;

	}

	friend std::ostream& operator<<(std::ostream& out, const StupidFunc & f ) {
		out << F( 8, 3, f.x_ ) << F( 8, 3, f.y_ ) << F( 8, 3, f.z_ ) << std::endl;
		return out;
	}

	core::Real x_;
	core::Real y_;
	core::Real z_;
};
