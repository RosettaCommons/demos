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

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Options.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PackMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <protocols/loops/ccd_closure.hh>
#include <protocols/relax_protocols.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>

#include <core/util/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
//#include <core/chemical/residue_io.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>

#include <core/mm/mm_torsion_library.hh>
#include <core/mm/mm_torsion_library.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/options/option.hh>
#include <core/options/after_opts.hh>

#include <core/util/basic.hh>

#include <core/io/database/open.hh>

#include <core/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//silly using/typedef


#include <core/util/Tracer.hh>

using core::util::T;
using core::util::Error;
using core::util::Warning;



using namespace core;
using namespace protocols;

using utility::vector1;

using io::pdb::dump_pdb;

void relax_test() {
 using namespace protocols::moves;
 using namespace scoring;

 pose::Pose pose;
 std::cerr << "READING 1rdg.pdb" << std::endl;
 io::pdb::pose_from_pdb( pose, "input/1rdg.pdb" ); // default is standard fullatom residue_set

	std::cerr << "SETUP SCORE FUNCTION" << std::endl;
 ScoreFunctionOP scorefxn( new ScoreFunction() );

 // aiming for standard packer weights
 scorefxn->set_weight( fa_atr, 0.80 );
 scorefxn->set_weight( fa_rep, 0.44 );
 scorefxn->set_weight( fa_sol, 0.65 );
 scorefxn->set_weight( fa_pair, 0.49 );
 scorefxn->set_weight( fa_dun, 0.56 );
 scorefxn->set_weight( rama, 0.2 );
 scorefxn->set_weight( hbond_lr_bb, 1.17 );
 scorefxn->set_weight( hbond_sr_bb, 1.17 );
 scorefxn->set_weight( hbond_bb_sc, 1.17 );
 scorefxn->set_weight( hbond_sc   , 1.10 );

	std::cerr << "SETING UP RELAX" << std::endl;
 protocols::ClassicRelax myrelax( scorefxn, pose );

	std::cerr << "RUNNING RELAX" << std::endl;
 myrelax.run();

	std::cerr << "DONE" << std::endl;
}

// void cpp_init() {
//   int argc = 3;
//   char ** argv = new char*[3];
//   argv[0] = new char[999];
//   argv[1] = new char[999];
//   argv[2] = new char[999];
//
//   std::strcpy(argv[0],"dummy");
//   std::strcpy(argv[1],"-database");
//   std::strcpy(argv[2],"/Users/sheffler/svn/branches/minirosetta_database");
//
//   std::cout << "MAKING ARGV:" << std::endl;
//   std::cout << "argv[0]: '" << argv[0] << "'" << std::endl;
//   std::cout << "argv[1]: '" << argv[1] << "'" << std::endl;
//   std::cout << "argv[2]: '" << argv[2] << "'" << std::endl;
//
//   core::init( argc, argv );
//
//  numeric::random::RandomGenerator::initializeRandomGenerators(
//     111111, numeric::random::_RND_TestRun_, "ran3");
//
// }

int main( int argc , char * argv [] ) {
	argc = 0;  argv = NULL;

	core::init("/Users/sheffler/svn/branches/minirosetta_database");
	numeric::random::RandomGenerator::initializeRandomGenerators(
		 111111, numeric::random::_RND_TestRun_, "ran3");

	relax_test();

	exit(0);


}

#define PYTHON
#ifdef  PYTHON

#include "boost/python.hpp"
#include "boost/python/suite/indexing/indexing_suite.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "boost/python/suite/indexing/map_indexing_suite.hpp"

pack::task::PackerTaskOP create_packer_task_as_reference( core::pose::Pose const & pose ) {
	using namespace core::pack::task;
	PackerTaskOP t = core::pack::task::TaskFactory::create_packer_task(pose);
	 //#ifdef SHEFFLERDEBUG
	 std::cerr << "relax.cc:198 (" << ")" << std::endl;
	 //#endif
	 t->initialize_from_command_line();//.restrict_to_repacking();
	 //#ifdef SHEFFLERDEBUG
	 std::cerr << "relax.cc:206 (" << ")" << std::endl;
	 //#endif

	 return t;
}

#include "core/pack/task/PackerTask_.hh"

using namespace utility::pointer;

struct Base : ReferenceCount {
	Base() {}
	virtual void foo() {std::cerr << "FOO!" << std::endl;}
};

struct Derived : Base {
	Derived() : Base() {}
	void foo() {std::cerr << "DERIVED FOO!" << std::endl;}
};

void base_func( owning_ptr<Base> b) {
	std::cerr << "GOT BASE" << std::endl;
	b->foo();
}

void derived_func( owning_ptr<Derived> d) {
	std::cerr << "GOT DERIVED" << std::endl;
	d->foo();
}

BOOST_PYTHON_MODULE_INIT(_test) {
	namespace bp = boost::python;
	using namespace core::pack::task;

	bp::class_<Base   ,owning_ptr<Base>    >("Base"   );
	bp::class_<Derived,owning_ptr<Derived> >("Derived");
	bp::def("base_func", &base_func);
	bp::def("derived_func", &derived_func);
	bp::implicitly_convertible< owning_ptr<Derived>, owning_ptr<Base> >();

//  bp::class_< PackerTask, PackerTaskOP, boost::noncopyable >("PackerTask",bp::no_init);

//  bp::class_< PackerTask_ >("PackerTask_");

	bp::def( "create_packer_task", &core::pack::task::TaskFactory::create_packer_task  );


	// bp::def( "relax_test"     , &relax_test       );
	//
	// bp::class_<     std::vector <std::string> > ("VecStr" ).def(bp::vector_indexing_suite< std::vector<std::string> >());
	// bp::class_< utility::vector1<std::string> > ("Vec1Str").def(bp::vector_indexing_suite< utility::vector1<std::string> >());
	// bp::class_<     std::vector <int> > ("Vec1Int").def(bp::vector_indexing_suite< std::vector<int> >() );
	// bp::class_< utility::vector1<int> > ("VecInt" ).def(bp::vector_indexing_suite< utility::vector1<int> >() );



}

#endif


