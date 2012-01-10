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


using namespace core;
using namespace core::optimization;

///////////////////////////////////////////////////////////////////////////////

class Histogram
{

	// todo: figure out how to zero a utility::vector1 after adding new values to
	// it. deal with accessing variables outside of the min/max range, especially
	// negative values.
public:
	Histogram(
		Real bin_size
	) : bin_size_( bin_size ) {
		max_val_ = 0;
		min_val_ = 0;
		npoints_ = 0;
	}

	void add_value( Real value ) {
		if ( value > max_val_ ) {
			densities_.resize( get_bin_index_( value ) );
			max_val_ = value;
		}
		int bin_index = get_bin_index_( value );
		++densities_[ bin_index ];
		++npoints_;
	}

	inline Real get_count( int index ) const {
		// return (float) (densities_[ index ] / npoints_);
		return densities_[ index ];
	}

	inline Real get_density( int index ) const {
		// return (float) ( this->get_density(index) / this->npoints() );
		return (float) this->get_count(index) / this->npoints();
	}

	Real bin_size() const {
		return bin_size_;
	}

	int nbins() const {
		return densities_.size();
	}

	int npoints() const {
		return npoints_;
	}

private:

	inline int get_bin_index_( Real value ) {
		// stupid indexing routine
		int i;
		Real val = bin_size_;
		for ( i = 1; val < value; ++i ) {
			val += bin_size_;
		}

		return i;
	}

	int  npoints_;
	Real min_val_;
	Real max_val_;
	Real bin_size_;

	utility::vector1< unsigned int > densities_;

	friend std::ostream& operator<<(std::ostream& out, const Histogram & h ) {
		out << A( 8, "bin_idx" )
				<< A( 8, "value"   )
				// << A( 8, "count"   )
				// << A( 8, "npoints" )
				<< A( 8, "density" )
				<< std::endl;

		for ( int i = 1; i <= h.nbins(); i++ ) {
			out << I( 8, i )
					<< F( 8, 3, i * h.bin_size() )
					// << F( 8, 3, h.get_count(i) )
					// << I( 8,    h.npoints()      )
					<< F( 8, 3, h.get_density(i) )
					<< std::endl;
		}
		return out;
	}

};

// mixture func to be fit to a histogram of densities.
class DensityFit : public core::optimization::Multifunc
{

public:
	DensityFit (
		Histogram h
	) : h_( h )
	{
		estimated_densities.resize( h_.nbins() );
	}

	inline
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const
	{
		// setup vars
		Real gaussian_param_  = vars[1];
		Real exp_param_       = vars[2];
		Real mixture_param_   = vars[3];
		Real anchor           = 0;
		Real sum_squared_diff = 0;

		// derivative with respect to exponential_param (b):
		// f = b * std::exp( -1 * b * x )
		// f'db = b * std::exp( -1 * b * x ) - x * std::exp( -1 * b * x )

		for ( int i = 1; i <= h_.nbins() ; ++i ) {
			Real x = i * h_.bin_size();

			Real const sqrt_2pi = 2.50662721600161f;
			Real r = std::abs( x - anchor );
			Real exp_score      =      mixture_param_  * std::exp( -1 * exp_param_ * r );
			Real gaussian_score = (1 - mixture_param_) *
																(1 / (gaussian_param_ * sqrt_2pi)) *
																std::exp( -0.5 * gaussian_param_ * r * r );
			Real density_calculated = exp_score + gaussian_score;
			sum_squared_diff += ( density_calculated - h_.get_density(i) ) *
													( density_calculated - h_.get_density(i) );

	}

		return std::sqrt( sum_squared_diff );
	}

	inline
	core::Real
	func( core::optimization::Multivec const & vars ) const
	{
		// setup vars
		Real gaussian_param_  = vars[1];
		Real exp_param_       = vars[2];
		Real mixture_param_   = vars[3];
		Real anchor           = 0;
		Real sum_squared_diff = 0;

		// derivative with respect to exponential_param (b):
		// f = b * std::exp( -1 * b * x )
		// f'db = b * std::exp( -1 * b * x ) - x * std::exp( -1 * b * x )


		for ( int i = 1; i <= h_.nbins() ; ++i ) {
			Real x = i * h_.bin_size();

			Real const sqrt_2pi = 2.50662721600161f;
			Real r = std::abs( x - anchor );
			Real exp_score      =      mixture_param_  * std::exp( -1 * exp_param_ * r );
			Real gaussian_score = (1 - mixture_param_) *
																(1 / (gaussian_param_ * sqrt_2pi)) *
																std::exp( -0.5 * gaussian_param_ * r * r );
			Real density_calculated = exp_score + gaussian_score;

			sum_squared_diff += ( density_calculated - h_.get_density(i) ) *
													( density_calculated - h_.get_density(i) );
		}

		return std::sqrt( sum_squared_diff );
	}



	inline
	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const
	{
		// setup vars
		Real gaussian_param_ = vars[1];
		Real exp_param_      = vars[2];
		Real mixture_param_  = vars[3];
		Real anchor = 0;

		Real delta = 1.90735e-06;  // stupid way of numerically estimating derivaties. very stupid.
		for ( Size ii = 1; ii <= vars.size(); ++ii ) {
			Multivec temp_vars = vars;
			temp_vars[ii]     += delta;

			Real x_diff  = delta;
			Real y_diff  = this->func( temp_vars ) - this->func( vars );
			// Real y_diff = 0.1;
			dE_dvars[ii] = y_diff / x_diff;
		}
	}

	Histogram histogram() {
		return h_;
	}

private:
	Histogram h_;
	utility::vector1< Real > estimated_densities;
};
