// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/r_frag_quality.cc
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange


#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <core/kinematics/MoveMap.hh>

#include <protocols/abinitio/util.hh>
#include <core/scoring/rms_util.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/pose/Pose.hh>
#include <core/init.hh>
#include <numeric/angle.functions.hh>
#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>


using core::util::T;
using core::util::Error;
using core::util::Warning;

//static core::util::Tracer TR("core.fragment.ConstantLengthFragments.cxxtest");
std::string const TR("core.fragment.ConstantLengthFragments.cxxtest");

using namespace core;
using namespace fragment;
using namespace pose;
using namespace options;
using namespace options::OptionKeys;

// hacky test code ---  should live somewhere else
void steal_constant_length_frag_set_from_pose ( pose::Pose const& pose, ConstantLengthFragSet& fragset ) {
	Size len = fragset.max_frag_length();
	Size const nbb ( 3 ); //steal phi, psi and omega
	kinematics::MoveMap dummy_movemap; //is ignored right now
	for ( Size pos = 1; pos <= pose.total_residue() - len + 1; ++pos ) {
		FragDataOP frag_raw = new FragData;
		for ( Size i = 1; i<= len; i++ ) {
			frag_raw->add_residue( new BBTorsionSRFD( nbb, pose.secstruct(pos), oneletter_code_from_aa(pose.residue( pos ).aa() ) ) );
		};
		FrameOP frame = new Frame( pos, len );
		frag_raw->steal( pose, *frame );
		frame->add_fragment ( frag_raw );
		fragset.add_frame( frame );
	};
}

Real compare_cartesian_rmsd( Pose const &orig_frag, Pose const &pred_frag ) {
	return scoring::rmsd_with_super( orig_frag, pred_frag, scoring::is_protein_backbone );
}

inline Real sqr ( Real x ) {
	return x*x;
}

Real compare_torsion_rmsd( Pose const &orig_frag, Pose const &pred_frag ) {
	Real err ( 0.0 );
	for ( Size pos = 1; pos <= orig_frag.total_residue(); ++pos ) {
		for ( Size dof = 1; dof <= 2; ++dof ) { //check phi and psi
			Real orig = orig_frag.torsion( id::TorsionID( pos, id::BB, dof ) );
			Real pred = pred_frag.torsion( id::TorsionID( pos, id::BB, dof ) );
			std::cout << orig << ' ' << pred << ' ' << orig-pred << " "
								<< numeric::nearest_angle(orig-pred,0.0)  << std::endl;
			err += sqr( numeric::nearest_angle(orig-pred,0.0) );
		}
	}
	return sqrt( err / 2*orig_frag.total_residue() );
}

Real compare_frags( Pose const &orig_frag, Pose const &pred_frag ) {
	return compare_cartesian_rmsd( orig_frag, pred_frag );
	//compare_torsion_rmsd( orig_frag, pred_frag );
}

int main( int argc, char** argv ) {

	core::init( argc, argv );

	kinematics::MoveMap move_all;
	int const frag_length = option[ abinitio::number_3mer_frags ];
	std::string const frag_file( option[  in::file::s ][1] );
	std::string const native_pdb ( option[ in::file::native ]() );
	ConstantLengthFragSet predicted_frags( frag_length );
	predicted_frags.read_fragment_file( frag_file );

	Pose native;
	//read it
	io::pdb::pose_from_pdb( native, native_pdb );
	ConstantLengthFragSet native_frags( frag_length );
	steal_constant_length_frag_set_from_pose( native , native_frags );

	Pose orig_frag; //an original fragment
	Pose pred_frag; //a predicted fragment

	FrameIterator fr_nat = native_frags.begin();

	//initialize poses
	fr_nat->fragment_as_pose( 1, orig_frag);
	fr_nat->fragment_as_pose( 1, pred_frag);

	for ( FrameIterator
		efr_nat=native_frags.end(),
		fr_pred=predicted_frags.begin(),
		efr_pred=predicted_frags.end();
	fr_nat!=efr_nat && fr_pred!=efr_pred;
	++fr_pred, ++fr_nat
			) {
		fr_nat->fragment( 1 ).apply( orig_frag, 1, frag_length );
		for ( Size i=1; i<=fr_pred->nr_frags(); i++ ) {
			fr_pred->fragment( i ).apply( pred_frag, 1, frag_length );
			std::cout << fr_pred->start() << i << compare_frags( orig_frag, pred_frag) << std::endl;
		}
	}
}
