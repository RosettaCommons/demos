// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file demo/phil_debug2.cc
/// @brief

/// @author James Thompson

#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>
#include <core/types.hh>

#include <core/init.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/Pose.hh>

#include <core/io/pdb/pose_io.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	// options, random initialization
	core::init( argc, argv );
	using namespace core::options;
	using namespace core::options::OptionKeys;

	// input information for debugging!
	std::string filename = "1ubi.pdb";

	core::pose::Pose pose;
	core::io::pdb::pose_from_pdb( pose, filename );
	// Phil: a debug build of this program should fail an assertion. If you move
	// the N-terminal Nitrogen away from (0,0,0), then it the assertion no longer
	// fails. Here's the stack trace from my macbook:

	// #0  0x9003d66c in kill ()
	// #1  0x9010e8cf in raise ()
	// #2  0x9010d422 in abort ()
	// #3  0x01e9e3d5 in __eprintf () at src/core/kinematics/Jump.hh:43
	// #4  0x01eaeb50 in numeric::xyzVector<double>::normalize (this=0xbfffe7a8) at src/numeric/xyzVector.hh:674
	// #5  0x052ce61b in core::kinematics::Stub::from_four_points (this=0xbfffe890, center=@0x126903d4, a=@0x126903d4, b=@0x12690694, c=@0x126906f4) at src/core/kinematics/Stub.cc:60
	// #6  0x054d8946 in core::kinematics::Stub::Stub (this=0xbfffe890, center=@0x126903d4, a=@0x126903d4, b=@0x12690694, c=@0x126906f4) at src/core/kinematics/Stub.hh:88
	// #7  0x052df483 in core::kinematics::Atom_::get_stub (this=0x126903c0) at src/core/kinematics/tree/Atom_.cc:414
	// #8  0x052e28c2 in core::kinematics::JumpAtom::update_internal_coords (this=0x126903c0, stub=@0x5584260, recursive=true) at src/core/kinematics/tree/JumpAtom.cc:168
	// #9  0x052c2ada in core::kinematics::AtomTree::update_internal_coords (this=0xbffff434) at src/core/kinematics/AtomTree.cc:1199
	// #10 0x052c50fc in core::kinematics::AtomTree::replace_tree (this=0xbffff434, new_root=0x126903c0, from_xyz=true) at src/core/kinematics/AtomTree.cc:165
	// #11 0x052f3c8b in core::conformation::Conformation::append_residue (this=0xbffff150, new_rsd_in=@0x126901f0, attach_by_jump=false, root_atomno=0, anchor_id=@0x5583c60) at src/core/conformation/Conformation.cc:398
	// #12 0x052f3d7c in core::conformation::Conformation::append_residue_by_bond (this=0xbffff150, new_rsd=@0x126901f0, build_ideal_geometry=false, residue_connection_index=0, anchor_pos=0, anchor_residue_connection_index=0) at src/core/conformation/Conformation.cc:276
	// #13 0x05381c06 in core::pose::Pose::append_residue_by_bond (this=0xbffff148, new_rsd=@0x126901f0, build_ideal_geometry=false, connection=0, anchor_residue=0, anchor_connection=0) at src/core/pose/Pose.hh:262
	// #14 0x0525c3c1 in core::io::pdb::FileData::build_pose_as_is (this=0xbffff03c, pose=@0xbffff148, residue_set=@0x1231a630) at src/core/io/pdb/file_data.cc:373
	// #15 0x0525c842 in core::io::pdb::FileData::build_pose (this=0xbffff03c, pose=@0xbffff148, residue_set=@0x1231a630) at src/core/io/pdb/file_data.cc:83
	// #16 0x05260fa0 in core::io::pdb::pose_from_pdb (pose=@0xbffff148, residue_set=@0x1231a630, filename=@0xbffff798) at src/core/io/pdb/pose_io.cc:95
	// #17 0x0526108a in core::io::pdb::pose_from_pdb (pose=@0xbffff148, filename=@0xbffff798) at src/core/io/pdb/pose_io.cc:108
	// #18 0x000027af in main (argc=3, argv=0xbffff808) at demo/james/phil_debug1.cc:51

	pose.dump_pdb("debug1.pdb");
	return 0;
}
