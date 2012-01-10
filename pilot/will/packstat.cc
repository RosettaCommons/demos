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
//#include <core/options/option.hh>

#include "devel/packstat/types.hh"
#include "devel/packstat/SimplePDB_Atom.hh"
#include "devel/packstat/SimplePDB.hh"
#include "devel/packstat/io.hh"
#include "devel/packstat/AtomRadiusMap.hh"
#include "devel/packstat/compute_sasa.hh"
#include "devel/packstat/sasa_dot_locations.hh"

#include "core/init.hh"
#include "core/types.hh"
#include "core/id/AtomID_Map.hh"
#include <core/options/option.hh>
#include <core/chemical/AtomType.hh>

#include "core/pose/Pose.hh"
#include "core/io/pdb/pose_io.hh"
#include "core/scoring/sasa.hh"

#include "core/util/Tracer.hh"

#include "utility/vector1.hh"
#include "utility/file/FileName.hh"
#include "utility/io/izstream.hh"
#include "utility/io/ozstream.hh"

#include "numeric/xyz.functions.hh"
#include "numeric/random/random.hh"

#include "ObjexxFCL/formatted.io.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctype.h>

core::util::Tracer TRps("packstat");

using namespace devel::packstat;

void test_io()
{
  using namespace devel::packstat;
  using namespace std;

  istringstream iss("ATOM     12  N   GLU A   2     -13.565  31.875  -5.182  1.00 51.33           N\n");
  SimplePDB_Atom atom;
  iss >> atom;
  cout << atom << endl;

  SimplePDB pdb;
  ifstream in("test_in.pdb");
  in >> pdb;

  cout << pdb;


}

void test_spheres(std::string fname)
{
  using namespace devel::packstat;
  using namespace std;

  AtomRadiusMap arm;
  SimplePDB pdb;
  utility::io::izstream in(fname.c_str());
  in >> pdb;

	// for( SPAtomIter i = pdb.atoms().begin(); i != pdb.atoms().end(); ++i ) {
	// 	TRps << *i << " " << arm.get_radius(i->type,i->res) << std::endl;
	// }

	// TRps << fname << std::endl;
  Spheres spheres( pdb.get_spheres(arm) );
  for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
     TRps << *i << std::endl;
  }

}


void test_sasa(std::string fname)
{
  using namespace devel::packstat;
  using namespace std;
	using namespace core;

	pose::Pose pose;
	io::pdb::pose_from_pdb(pose,fname);

	utility::vector1< Real > radii;
	{
		chemical::AtomTypeSet const & atom_set( pose.residue(1).atom_type_set() );
		Size const SASA_RADIUS_INDEX( atom_set.extra_parameter_index( "SASA_RADIUS" ) );
		radii.resize( atom_set.n_atomtypes() );
		for ( Size i=1; i<= radii.size(); ++i ) {
			chemical::AtomType const & at( atom_set[i] );
			radii[i] = atom_set[i].extra_parameter( SASA_RADIUS_INDEX );
		}
	}

	Spheres pose_spheres;
	{
		id::AtomID_Map<Real> atom_sasa;
		utility::vector1<Real> rsd_sasa;
		scoring::calc_per_atom_sasa(pose,atom_sasa,rsd_sasa,3.0);
		for( size_t i = 1; i <= pose.total_residue(); ++i ) {
			for( size_t j = 1; j <= pose.residue(i).natoms(); ++j ) {
				TRps << "3.0 " << i << " " << j << " " << atom_sasa[i][j] << std::endl;
				conformation::Atom const & a( pose.residue(i).atom(j) );
				pose_spheres.push_back( Sphere( a.xyz(), radii[a.type()] ) );
			}
		}
		TRps << std::endl;
	}
	{
		id::AtomID_Map<Real> atom_sasa;
		utility::vector1<Real> rsd_sasa;
		scoring::calc_per_atom_sasa(pose,atom_sasa,rsd_sasa,1.4);
		for( size_t i = 1; i <= pose.total_residue(); ++i ) {
			for( size_t j = 1; j <= pose.residue(i).natoms(); ++j ) {
				TRps << "1.4 " << i << " " << j << " " << atom_sasa[i][j] << std::endl;
			}
		}
		TRps << std::endl;
	}

	AtomRadiusMap arm;
  SimplePDB pdb;
  utility::io::izstream in(fname.c_str());
  in >> pdb;

	// TRps << fname << std::endl;
  Spheres spheres( pdb.get_spheres(arm) );
  // for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
  //   TRps << *i << std::endl;
  // }
	spheres = pose_spheres;

	SasaOptions opts;
	opts.probe_radii.push_back(3.0);
	opts.probe_radii.push_back(1.4);
	SasaResultOP sr = compute_sasa( spheres, opts );
	for( size_t p = 1; p <= 2; ++p ) {
		for( size_t i = 1; i <= spheres.size(); ++i ) {
			TRps << opts.probe_radii[p] << " " << i << " " << sr->sphere_sasa(i,p) << std::endl;
		}
	}

}

// not very interesting
// void test_cav_ball_overlap( devel::packstat::CavBalls & cbs ) {
// 	using namespace devel::packstat;
// 	for( PackstatReal delta = 0; delta < 1; delta += 0.1 ) {
// 		int count = 0;
// 		for( CavBallIter i = cbs.begin(); i != cbs.end(); ++i ) {
// 			CavityBall & cb1( *i );
// 			for( CavBallIter j = cbs.begin(); j != cbs.end(); ++j ) {
// 				CavityBall & cb2( *j );
// 				PackstatReal dist = distance( cb1.xyz(), cb2.xyz() );
// 				if( 0 == dist ) continue;
// 				if( dist + cb1.radius() <= cb2.radius()+delta ) {
// 					count++;
// 				}
// 			}
// 		}
// 		TRps << delta << " " << count << std::endl;
// 	}
// }

inline std::string base_name(const std::string& str)
{
  size_t begin = 0;
  size_t end = str.length();

  for (size_t i=0; i<str.length(); ++i) {
    if (str[i] == '/') begin = i+1;
  }

  return str.substr(begin,end);
}

std::string get_out_tag(std::string fname) {
	std::string base = base_name(fname);
	std::transform( base.begin(), base.end(), base.begin(), tolower );
	system( ("mkdir -p out/" + base.substr(1,2)).c_str() );
	std::string OUT_TAG = "out/" + base.substr(1,2) + "/" + base;
	return OUT_TAG;
}

struct OrderSphereOnX {
  bool operator()( devel::packstat::Sphere const & a, devel::packstat::Sphere const & b ) {
		return a.xyz.x() < b.xyz.x();
  }
};


void test_cav_balls( std::string fname ) {
	using namespace devel::packstat;
  using namespace std;
	using namespace core;
	using namespace ObjexxFCL::fmt;
	using namespace numeric;
	using namespace utility;

	string OUT_TAG = get_out_tag(fname);
	TRps << "OUT_TAG " << OUT_TAG << std::endl;

	AtomRadiusMap arm;
  SimplePDB pdb;
  utility::io::izstream in(fname.c_str());
  in >> pdb;

	TRps << fname << std::endl;
  Spheres spheres( pdb.get_spheres(arm) );
	std::sort( spheres.begin(), spheres.end(), OrderSphereOnX() );

	vector1< xyzVector<PackstatReal> > centers( pdb.get_res_centers() );

	SasaOptions opts;
	opts.prune_max_iters = 999;
	opts.prune_max_delta = 0;
	opts.num_cav_ball_layers = 10;
	opts.frac_cav_ball_required_exposed = 0.00;
	opts.area_cav_ball_required_exposed = 0.00;
	opts.surrounding_sasa_smoothing_window = 3;
	for( PackstatReal pr = 3.0; pr >= 0.4; pr -= 0.0333 ) opts.probe_radii.push_back(pr);
	opts.prune_cavity_burial_probe_radii.push_back( 1.6 );
	for( PackstatReal pr = 1.6; pr >= 0.1; pr -= 0.1    ) opts.prune_cavity_burial_probe_radii.push_back(pr);

	TRps << "compute MSAs" << std::endl;
	SasaResultOP sr = compute_sasa( spheres, opts );

	////////////////////////////////////////////////////////////////////////////////////////////////
	CavBalls cavballs = sr->cavballs;
	TRps << "pruning hidden cav balls " << cavballs.size() << std::endl;
	cavballs = prune_hidden_cavity_balls( cavballs, opts );
	TRps << "pruning exposed cav balls " << cavballs.size() << std::endl;
	cavballs = prune_cavity_balls( spheres, cavballs, opts );
	TRps << "compute cav ball volumes	" << cavballs.size() << std::endl;
	compute_cav_ball_volumes( cavballs, opts );
	TRps << "selecting representatives" << std::endl;
	CavBalls sel_cbs( select_cav_balls(cavballs,4.0) );

	ostringstream cav_info;
	for( CavBallIter i = cavballs.begin(); i != cavballs.end(); ++i ) {
		// PackstatReal maxa = i->radius() * i->radius() * 4.0 * 3.14159 * 1.00;
		// PackstatReal maxv = i->radius() * i->radius() * i->radius() * 4.0/3.0 * 3.14159 * 1.00;
		cav_info << F( 7, 3, i->radius() ) << " "
		         << F( 7, 3, i->exposed_radius ) << " "
		         << F( 7, 3, i->area     ) << " "
		         << F( 7, 3, i->vol      ) << " "
		         << std::endl;
		// if( i->area > maxa * 1.01 ) TRps << "area too high " << i->area << " " << maxa << std::endl;
		// if( i->vol  > maxv * 1.01 ) TRps << "vol. too high " << i->vol  << " " << maxv << std::endl;
	}
	utility::io::ozstream cav_info_file(((OUT_TAG+".cavities").c_str()));
	cav_info_file << cav_info.str();
	cav_info_file.close();
	cav_info.clear();


	////////////////////////////////////////////////////////////////////////////////////////////////
	// TRps << "SASA" << std::endl;
	// for( size_t i = 1; i <= spheres.size(); ++i ) {
	// 	TRps << I(5,i) << " ";
	// 	for( size_t p = opts.probe_radii.size(); p >= 1; --p ) {
	// 		TRps << F( 5, 2, sr->sphere_sasa(i,p) ) << " ";
	// 	}
	// 	TRps << std::endl;
	// }

	///////////////////////////////////////////////////////////////////////////////////////////////
	TRps << "writting stupid pdb to "+OUT_TAG+".pdb" << std::endl;
	ostringstream out;
	  for( SphereIter i = spheres.begin(); i != spheres.end(); ++i ) {
		int rnum = 0, anum = 0;
		PackstatReal occ = 0.0f;
	    out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
		+ I( 4, rnum ) + "    "
		+ F( 8, 3, i->xyz.x() ) + F( 8, 3, i->xyz.y() ) + F( 8, 3, i->xyz.z() )
		+ F( 6, 2, occ ) + ' ' + F( 5, 2, i->radius ) << std::endl;
	  }
	for( size_t cb = 1; cb <= cavballs.size(); ++cb ) {
		out << cavballs[cb].hetero_atom_line() << std::endl;
	}
	// for( size_t cb = 1; cb <= sel_cbs.size(); ++cb) {
	// 	out << sel_cbs[cb].hetero_atom_line(7) << std::endl;
	// }
	utility::io::ozstream outz((OUT_TAG+".pdb").c_str());
	outz << out.str();
	outz.close();
	out.clear();

	//////////////////////////////////////////////////////////////////////////////////////////////////
	TRps << "output_res_surrounding_sasa " << centers.size() << std::endl;
	ostringstream outss2;
	for( size_t i = 1; i <= centers.size(); ++i ) {
		outss2 << i << " ";
		output_surrounding_sasa( outss2, centers[i], spheres, sr, opts );
		outss2 << std::endl;
	}
	utility::io::ozstream outss2z((OUT_TAG+".res_sur_sasa").c_str());
	outss2z << outss2.str();
	outss2z.close();
	outss2.clear();

	TRps << "output_surrounding_sasa " << sel_cbs.size() << std::endl;
	ostringstream outss;
	for( size_t i = 1; i <= sel_cbs.size(); ++i ) {
		CavityBall & cb( sel_cbs[i] );
		outss << cb.radius_ << " "
					<< cb.exposed_radius << " "
					<< cb.area << " "
					<< cb.vol << " ";
		output_surrounding_sasa( outss, sel_cbs[i].xyz(), spheres, sr, opts );
		outss << std::endl;
	}
	utility::io::ozstream outssz((OUT_TAG+".cav_sur_sasa").c_str());
	outssz << outss.str();
	outssz.close();
	outss.clear();
}

numeric::xyzMatrix<PackstatReal> rand_rot() {
	using namespace numeric;
	using namespace numeric::random;
	xyzVector<PackstatReal> axis(uniform(),uniform(),uniform());
	while( axis.length() > 1 ) axis = xyzVector<PackstatReal>(uniform(),uniform(),uniform());
	return rotation_matrix<PackstatReal>( axis, uniform() * 2 * 3.14159 );
}

void test_sasa_dots() {
	using namespace utility;
	using namespace numeric;
	using namespace devel::packstat;
	using namespace std;

	vector1<xyzVector<PackstatReal> > sasa_dots = get_sasa_dot_locations();

	{
		ofstream out("sasa_dots_randomized.pdb");
		for( PackstatReal pr = 350.0; pr <= 600.0; pr += 10.0 ) {
			xyzMatrix<PackstatReal> rot = rand_rot();
			for( size_t i = 1; i <= sasa_dots.size(); ++i ) {
				xyzVector<PackstatReal> dot = rot * sasa_dots[i];
				int rnum = 0, anum = 0;
				PackstatReal occ = 0.0f;
			    out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
				+ I( 4, rnum ) + "    "
				+ F( 8, 3, dot.x()*pr ) + F( 8, 3, dot.y()*pr ) + F( 8, 3, dot.z()*pr )
				+ F( 6, 2, occ ) + ' ' + F( 5, 2, 5.0 ) << std::endl;
			}
		}
		out.close();
	}
	{
		ofstream out("sasa_dots.pdb");
		for( PackstatReal pr = 350.0; pr <= 600.0; pr += 10.0 ) {
			for( size_t i = 1; i <= sasa_dots.size(); ++i ) {
				xyzVector<PackstatReal> dot = sasa_dots[i];
				int rnum = 0, anum = 0;
				PackstatReal occ = 0.0f;
			    out << "ATOM  " + I( 5, ( anum ) ) + "  V   PRT Z"
				+ I( 4, rnum ) + "    "
				+ F( 8, 3, dot.x()*pr ) + F( 8, 3, dot.y()*pr ) + F( 8, 3, dot.z()*pr )
				+ F( 6, 2, occ ) + ' ' + F( 5, 2, 5.0 ) << std::endl;
			}
		}
		out.close();
	}

}

int main (int argc, char *argv[])
{

	core::init( argc, argv );

  using namespace devel::packstat;
  using namespace core::options;
  using namespace core::options::OptionKeys;
  using namespace utility;

  // test_io();

	// test_sasa_dots();

	if( option[ in::file::s ].user() ) {
  	vector1<file::FileName> files( option[ in::file::s ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
    	test_cav_balls( files[i] );
  	}
	} else if( option[ in::file::l ].user() ) {
  	vector1<file::FileName> files( option[ in::file::l ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
			utility::io::izstream list( files[i] );
			std::string fname;
			while( list >> fname ) {
				std::cerr << "'" << fname << "'" << std::endl;
    		test_cav_balls( fname );
			}
  	}
	}
	return 0;

}
