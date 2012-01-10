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

#include <core/types.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/VariantType.hh>

#include "core/scoring/ScoreFunctionFactory.hh"

#include <core/init.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>



#include "utility/pointer/owning_ptr.hh"
#include <utility/vector1.hh>

#include "ObjexxFCL/formatted.o.hh"


// C++ headers
#include <iostream>
#include <string>
#include <vector>





int main (int argc, char *argv[])
{

	core::init( argc, argv );
  using namespace core::scoring;

  ScoreFunctionOP sf( ScoreFunctionFactory::create_score_function(STANDARD_WTS) );


	return 0;

}
