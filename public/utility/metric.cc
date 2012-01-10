// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/utility/true_false.cc
/// @brief  demo of True/False query class which handles basic boolean logic
/// @author Will Sheffler (willsheffler@gmail.com)
/// @date   Thu Oct 19 19:49:23 2007
///

#include <cstdlib>
#include <iostream>
#include <vector>
#include <map>

#include "core/init.hh"
#include "core/types.hh"

#include "numeric/random/random.hh"

#include "utility/query/types.hh"
#include "utility/query/Metric.hh"

#include "utility/pointer/owning_ptr.hh"

#include "ObjexxFCL/string.functions.hh"

numeric::random::RandomGenerator RG(897987987);

using core::Real;
using namespace utility::query;
using std::string;
using std::map;
using std::vector;

void print_vals(Filter<Real> f) {
	std::cout << f.description() << std::endl;
	for(Real r = -5; r <= 5; r++) {
		if(f(r)) std::cerr << r << " True " << std::endl;
		else     std::cerr << r << " False" << std::endl;
	}
}

void print_vals(Metric<Real> f) {
	std::cout << f.description() << std::endl;
	for(Real r = -5; r <= 5; r++)
		std::cerr << r << " " << f(r) << std::endl;

}

void test_real() {
	std::cerr << "============== Real tests =============" << std::endl;

	Filter<Real> f;
	Metric0Param<Real> x(new ImplicitConverter<Real,Real>);

	print_vals( ~ (x*x*x < 0) );
	print_vals( (-1 + x*x + 1) == (x^2) ); // should be true!
	Metric<Real> m = ( (x + 3)*x - 2*x + 4/x );
	print_vals( m );

}



// some dummy structure link something from a silent file
struct SilentInfo {
	SilentInfo() {
		scores_["score"]  = RG.uniform() * 100.0;
		scores_["fa_atr"] = RG.uniform() * 100.0;
		scores_["fa_rep"] = RG.uniform() * 100.0;
	}
	Real get_score(std::string s) { return scores_[s]; }
	std::map<string,Real> scores_;
};


// must convert SilentInfo*,string pair to 'Real'
struct SilentInfoGetter : public Converter1Param<SilentInfo*,string,Real> {
	Real convert(SilentInfo * si, string s) { return si->get_score(s); }
	string description() { return "silent info score"; }
};



void test_fake_silent() {
	std::cout << "========================== si test ===============" << std::endl;

	// make a vector of silent data things
	std::vector<SilentInfo*> v(100);
	for( int ii = 0; ii <= 99; ++ii ) { v[ii] = new SilentInfo; }

	// these could be stuffed away in some namespace in an actual application
	// the user could just 'using namespace silent::query' and get the appropriate stuff
	SilentInfoGetter * sig = new SilentInfoGetter;
	Metric1Param<SilentInfo*,string> fa_atr( sig, "fa_atr" );
	Metric1Param<SilentInfo*,string> fa_rep( sig, "fa_rep" );
	Metric1Param<SilentInfo*,string> score ( sig, "score" );

	// make whatever compound query you want. supports +-*/ and ^ for exponent
	Metric<SilentInfo*> myscore = ( fa_atr + fa_rep/12 + score^2 ) ^ 0.3;
	std::cout << myscore.description() << std::endl;

	// putting an inequality ==, !=, <, <=, >, or >=automatically turns the
	// 'Metric' into a 'Filter' which can be combined w/ other filters with
	// &, | and ~ for not
	Filter<SilentInfo*> filt = myscore > 19.0 & fa_atr > 7;
	std::cout << filt.description() << std::endl;

	// a test...
	for( int ii = 0; ii <= 99; ++ii ) {
		SilentInfo * si( v[ii] );
		if( !filt(si) ) continue; // calling 'filt' returns a bool
		std::cerr << ii << " " << myscore(si) << std::endl; // calling myscore returns a Real
	}

}



int main (int argc, char * argv[])
{
	core::init(argc,argv);

	test_real();

	test_fake_silent();

	std::cerr << "DONE!" << std::endl;

	return 0;

}


