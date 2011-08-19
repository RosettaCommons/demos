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
/// @author will sheffler

#include <core/types.hh>
#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <utility/exit.hh>

#include "core/io/serialization/serialize_pose.hh"
#include "core/pose/Pose.hh"
#include "core/io/pdb/pose_io.hh"


std::string DEFAULT_DB = "minirosetta_database_rev19451.zip";
using namespace core::io::serialization;
using std::string;

void test_serialization()
{
	const size_t N = 99999;
	char buf_b[N];
	BUFFER b((char*)(&buf_b),N);
	char c = '_';
	unsigned int i = -1;
	string s;
	float f;
	double d;
	char buf[120];
	string name = "DUMMY";

	write_binary('C',b);
	read_binary ( c ,b);
	std::cerr << "'" << c << "' should be C" << std::endl;

	write_binary(187u,b);
	read_binary ( i ,b);
	std::cerr << "'" << i << "' should be 187" << std::endl;

	write_binary(string("fubar!! bla"),b);
	write_binary(37.53f,b);
	write_binary(43234.8453,b);

	read_binary (s,b);
	std::cerr << "'" << s << "' should be 'fubar!! bla'" << std::endl;

	read_binary (f,b);
	std::cerr << "'" << f << "' should be '37.53f'" << std::endl;

	read_binary (d,b);
	std::cerr << "'" << d << "' should be '43234.8453'" << std::endl;



	core::pose::Pose p,p2;
	core::io::pdb::pose_from_pdb(p,"/Users/sheffler/svn/branches/mini/demo/phil/input/test_in.pdb");
	core::io::pdb::pose_from_pdb(p2,"/Users/sheffler/svn/branches/mini/demo/phil/input/test_in.pdb");
	p2.set_phi(7,-23.456);
	std::cerr << p.phi(7) << " " << p2.phi(7) << std::endl;

	std::cerr << "serializing pose!!" << std::endl;
	write_binary(p,b);

	std::cerr << "unserializing pose!!" << std::endl;
	read_binary(p2,b);

	std::cerr << p.phi(7) << " " << p2.phi(7) << std::endl;

}


int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );

	test_serialization();

}

