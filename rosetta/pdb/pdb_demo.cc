// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/io/pdb/file_data.cc
///
/// @brief
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
//#include <core/chemical/residue_io.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>



#include <core/io/database/open.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/file_data.hh>

#include <core/util/OStream.hh>
#include <core/util/Tracer.hh>
#include <core/util/basic.hh>
#include <core/init.hh>
#include <core/types.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>


//silly using/typedef

using namespace core;
using namespace core::util;

using utility::vector1;
using core::util::T;
using core::util::Error;
using core::util::Warning;


typedef std::map< std::string, std::map< std::string, numeric::xyzVector< Real > > > Coords;

typedef vector1< std::string > Strings;

///////////////////////////////////////////////////////////////////////////////
// some silly helper routines:
//
///////////////////////////////////////////////////////////////////////////////

std::string readFile(std::string fname)
{
	Size fsize;
	std::string res;

	std::ifstream file(fname.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
	if( file.is_open() )  {
		fsize = file.tellg();
		res.resize( fsize );
		file.seekg(0, std::ios::beg);
		file.read(&res[0], fsize);
		file.close();
	}
	else std::cout << "file not found!";
	return res;
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] )
{

	using namespace core;

	core::init(argc, argv);

	chemical::ResidueTypeSetCAP residue_set
		( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

	T("pdb_demo") << "---------------------\n";

	T("pdb_demo") << "\n\n\n\n\n>>>>> Reading file...\n";
	//std::string data = readFile("test_in.pdb");
	std::string data = readFile("test_in.pdb");

	T("pdb_demo") << ">>>>> Creating file data...\n";
	core::io::pdb::FileData fd = core::io::pdb::PDB_DReader::createFileData(data);

	T("pdb_demo") << ">>>>> Building Pose...\n";

	//T("pdb_demo") << fd;

		core::pose::Pose pose;
	fd.build_pose(pose, *residue_set);
	T("pdb_demo") << " pose.total_residue()=" << pose.total_residue() << "\n";
	//T("pdb_demo") << fd;
	T("pdb_demo") << "Back to PDB now... -----------------\n";

	core::io::pdb::FileData fd1;
	fd1.init_from_pose(pose);
	//std::cout << io::pdb::PDB_DReader::createPDBData(fd1);

	pose.dump_pdb("output.pdb");


	using namespace options::OptionKeys;

	std::vector<int> A;  A.push_back(1);  A.push_back(2);  A.push_back(3);  A.push_back(5);
	utility::vector1<int> B;  B.push_back(10);  B.push_back(20);  B.push_back(30);  B.push_back(45);
	std::map<int, std::string> M;  M[1]="one";  M[2]="two";  M[3]="1+2";

	T("Demo") << "vector:" << A << " vector1:" << B << " map:" << M << "\n";

	T("Error", 10) << "Some error here!!!\n";
	T("core.pose") << "Some core pose message\n";
	T("core") << "Some core message\n";

	Error() << "Some error test...\n";
	Warning() << "Some warning test...\n";

	std::cout << "Done! -------------------------------\n";
}

