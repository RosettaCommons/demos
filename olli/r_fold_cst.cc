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


#include <core/types.hh>
#include <core/init.hh>

#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>
#include <core/options/after_opts.hh>
#include <core/util/Tracer.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentFileData.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>


#include <protocols/relax_protocols.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>

#include <utility/vector1.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <protocols/viewer/viewers.hh>

core::util::ATracer("r_fold_cst") tr;

using namespace core;
using namespace protocols;
using namespace fragment;
using namespace abinitio;

#include "helper_code_r_trjconv.cc"

using namespace core;

#include <devel/simple_options/option.hh>

using namespace devel::option;

class ThisApplication {
public:
	ThisApplication();
	static void register_options();
	void process_decoy( pose::Pose &pose, std::string tag );
	void add_constraints( pose::Pose &pose );
	void do_rerun();
	void fold();
	void run();
	void setup();
private:
	std::ofstream score_file_;
	pose::PoseOP rmsd_pose_;
	pose::Pose native_pose_;
	PCA_OP pca_;
	scoring::ScoreFunctionOP scorefxn_;
	scoring::EMapVector output_weights_;
	std::string sequence_;
	scoring::constraints::ConstraintSetOP cstset_;
};

ThisApplication::ThisApplication() :
	rmsd_pose_( NULL),
	pca_( NULL ),
	cstset_( NULL )
{}


using namespace options;
using namespace options::OptionKeys;


#define OPT(akey)																					\
	core::options::option.add_relevant( akey )

#define NEW_OPT(akey,help,adef)								      		  \
	core::options::option.add( akey , help ).def( adef ); 	\
	OPT( akey )

#define OPT_KEY( type, key )                                    \
	namespace core { namespace options { namespace OptionKeys {		\
				type##OptionKey const key( #key );											\
	} } }

OPT_KEY( Boolean, rerun )
OPT_KEY( Boolean, steal )
OPT_KEY( Boolean, start_extended )
OPT_KEY( File, pca )
OPT_KEY( File, rmsd_target )
OPT_KEY( File, sf )
OPT_KEY( Boolean, viol )
OPT_KEY( Integer, viol_level )
OPT_KEY( String, viol_type )
OPT_KEY( StringVector, tag_selector )

void ThisApplication::register_options() {
	using namespace core::options;
	using namespace core::options::OptionKeys;
	OPT( in::file::native );
	OPT( in::file::silent ); // input silent file
	OPT( out::file::silent );
	OPT( in::file::frag3 );
	OPT( in::file::frag9 );
	OPT( in::file::fasta );
	OPT( constraints::cst_file );
	//	OPT( constraints::cst_weight );
	OPT( out::nstruct );
	NEW_OPT( OptionKeys::rerun, "go through intput structures and evaluate ( pca, rmsd, cst-energy )", false );
	NEW_OPT( steal, "use fragments of native (starting) structure", true);
	NEW_OPT( start_extended, "start from extended structure (instead of native)", false );
	NEW_OPT( pca, "compute PCA projections", "");
	NEW_OPT( rmsd_target, "compute rmsd against this structure", "" );
	NEW_OPT( sf, "filename for score output", "score.fsc" );
	NEW_OPT( viol, "show violations", true );
	NEW_OPT( viol_level, "how much detail for violtion output", 1 );
	NEW_OPT( viol_type, "work only on these types of constraints", "");
	NEW_OPT( tag_selector, "work only on these tag(s)","");
	FoldConstraints::register_options();
}

void ThisApplication::setup() {
	core::options::option[ options::OptionKeys::in::path::database ].def( "~/minirosetta_database");
	score_file_.open( std::string(option[ sf ]()).c_str() );
	if ( score_file_.fail() ) {
		trWarning << "unable to open score_file "+std::string( option[ sf ]() )+" for output\n";
	}

	//	native_pose_.clear();
	io::pdb::pose_from_pdb( native_pose_, option[ in::file::native ]() );
	switch_to_residue_type_set( native_pose_, chemical::CENTROID );

	if ( option [ in::file::fasta ].specified() ) {
		sequence_ = read_fasta ( option[ in::file::fasta ]() );
	} else {
		sequence_ = native_pose_.sequence();
	}

	if ( option[ rmsd_target ].user() ) {
		rmsd_pose_ = new pose::Pose;
		rmsd_pose_->clear();
		io::pdb::pose_from_pdb( *rmsd_pose_, option[ rmsd_target]() );
	}

	if ( option[ pca ].user() ) {
		pca_ = new PCA;
		pca_->read_eigvec_file( option[ pca ](), native_pose_, 2 );
		pca_->show( trTrace );
	}

	scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "score3" );
	scorefxn_->set_weight( scoring::atom_pair_constraint, option [ constraints::cst_weight ] );
	output_weights_ = scorefxn_->weights();
	if ( rmsd_pose_ ) output_weights_[ core::scoring::rms ] = 1.0;
	if ( pca_ ) output_weights_[ core::scoring::fa_atr ] = 1.0;
	if ( pca_ ) output_weights_[ core::scoring::fa_dun ] = 1.0;

}

void ThisApplication::process_decoy( pose::Pose &pose, std::string tag ) {
	//pose.dump_pdb("test.pdb");
	if ( option [ viol ] ) {
		pose.constraint_set()->show_violations(  std::cout, pose, option[ viol_level ] );
	};

	( *scorefxn_ )( pose );
	scorefxn_->show( trInfo, pose );

	if ( rmsd_pose_ ) { // calculate RMSD
		trInfo << "compute RMSD for " << tag << "\n";
		core::Real rmsd = core::scoring::CA_rmsd( *rmsd_pose_, pose );
		pose.energies().total_energies()[ core::scoring::rms ] = rmsd;
	}

	if ( pca_ ) { // compute projections to PCA vectors
		PCA::ProjectionVector proj;
		pca_->eval( pose, proj );
		trInfo << "PCA  " << proj[1] << " " << proj[2] << std::endl;
		pose.energies().total_energies()[ core::scoring::fa_atr ] = proj[1];
		pose.energies().total_energies()[ core::scoring::fa_dun ] = proj[2];
	}

	if ( score_file_.good() ) {
		pose.energies().total_energies().show_if_nonzero_weight( score_file_, output_weights_ );
		( score_file_ ) << " " << tag;
		( score_file_ ) << std::endl;
	}
}

void ThisApplication::add_constraints( pose::Pose& pose ){
	using namespace core::scoring::constraints;
	if ( ! cstset_ ) {
		cstset_ = ConstraintIO::read( option[ constraints::cst_file ], new ConstraintSet, pose );
	} else {
		pose.constraint_set( cstset_ );
	}
}

void ThisApplication::do_rerun() {
	using namespace core;
	using namespace io::silent;
	using namespace pose;

	//read silent file for input
	ProteinSilentFileData sfd;
	sfd.read_file( option [ in::file::silent ]() );
	// add native structure to the list
	sfd.add_structure( new ProteinSilentStruct( native_pose_, "NATIVE", false ) );  /// !!! VERY BAD IF NOT IDEALIZED !!! */
	// run thru all structures
	for ( ProteinSilentFileData::const_iterator it=sfd.begin_const(), eit=sfd.end_const(); it!=eit; ++it ) {
		Pose pose;
		std::string tag = it->decoy_tag();
		if ( option[ tag_selector ].user() == 0 || std::find( option[ tag_selector ]().begin(), option[ tag_selector ]().end(), tag ) != option[ tag_selector ]().end() ) {
			it->fill_pose( pose, *chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID ));
			if ( tag == "NATIVE" ) pose=native_pose_; // replace structure with NATIVE so that we don't suffer from non-idealized stuff
			add_constraints( pose );
			trInfo << tag << " " ;
			process_decoy( pose, tag );
		}
	}
}

void ThisApplication::fold() {
	// ==========================================================================
	///  --------- fold()-specific setup --------------------------
	// ==========================================================================
	ConstantLengthFragSetOP fragset3mer = new ConstantLengthFragSet( 3 );
	ConstantLengthFragSetOP fragset9mer = new ConstantLengthFragSet( 9 );
	fragset3mer->read_fragment_file( option [ in::file::frag3 ]);
	fragset9mer->read_fragment_file( option [ in::file::frag9 ], 25 );

	if ( option [ steal ] ) {
		steal_constant_length_frag_set_from_pose( native_pose_, *fragset3mer );
		steal_constant_length_frag_set_from_pose( native_pose_, *fragset9mer );
	};

	pose::Pose extended_pose;
	make_pose_from_sequence_( sequence_,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )),
		extended_pose
	);

	// make extended chain
	for ( Size pos = 1; pos <= extended_pose.total_residue(); pos++ ) {
		extended_pose.set_phi( pos, -45 );
		extended_pose.set_psi( pos, -45 );
		extended_pose.set_omega( pos, 180 );
	}

	// requires that the sequences match at the beginning -- > use sequence alignment later
	if ( ! option [ start_extended ] ) {
		trInfo << " *** use native structure as starting template *** \n";
		// determine length of segment to copy from native
		Size seg_len = std::min(extended_pose.total_residue(), native_pose_.total_residue() );
		fragment::Frame long_frame(1, seg_len);

		//create apropriate length FragData object
		FragData frag; //there should be some kind of factory to do this.
		Size nbb ( 3 ); //3 backbone torsions to steal
		for ( Size pos = 1; pos<= seg_len; pos++ ) {
			frag.add_residue( new BBTorsionSRFD( nbb, native_pose_.secstruct( pos ), 'X' ) );
		};
		// get torsion angles from native pose
		frag.steal(native_pose_, long_frame);

		// apply native torsions to extended structue
		frag.apply(extended_pose, long_frame);
	}; // if option [ start extended ]

	extended_pose.dump_pdb( "starting_structure.pdb", trInfo );
	add_constraints( extended_pose );



	// make a MoveMap
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	movemap->set_bb( true );



	FoldConstraints abinitio_protocol( fragset3mer, fragset9mer, movemap );
	abinitio_protocol.init( extended_pose );
	//  abinitio_protocol.set_constraint_weight( option[ constraints::cst_weight ]);

	// ==========================================================================
	// ==========================================================================
	// ==========================================================================

	//
	//     Start running abinitio
	//

	using protocols::jobdist::BasicJob;
	using protocols::jobdist::BasicJobOP;
	using protocols::jobdist::PlainSilentFileJobDistributor;
	utility::vector1< BasicJobOP > input_jobs;
	int const nstruct = std::max( 1, option [ out::nstruct ]() );
	BasicJobOP job = new BasicJob("classic_abinitio_relax", nstruct);
	input_jobs.push_back( job );
	PlainSilentFileJobDistributor< BasicJobOP > jobdist( input_jobs );
	BasicJobOP curr_job, prev_job;
	int curr_nstruct, num_structures_processed = 0;
	jobdist.startup();
	while ( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		std::cout << "Starting " << curr_job->output_tag(curr_nstruct) << " ...\n";

		pose::Pose fold_pose ( extended_pose );

		//    protocols::moves::MonteCarlo& mc = abinitio_protocol.mc();
		//    protocols::viewer::add_monte_carlo_silent_viewer( mc, "test_mc_out", false );
		int start_time = time(NULL);

		abinitio_protocol.apply( fold_pose );
		fold_pose.constraint_set()->show_violations(  trInfo, fold_pose, option [ viol_level ] );

		int end_time   = time(NULL);
		std::cout << "TIMEFORABINITIO: " << end_time - start_time << std::endl;

		process_decoy( fold_pose, curr_job->output_tag( curr_nstruct ) );

		tr.Debug << "RMSD in pose is " << fold_pose.energies().total_energies()[ scoring::rms ] << "\n";
		bool fullatom( false );
		jobdist.dump_pose( curr_nstruct, fullatom, fold_pose );
		prev_job = curr_job;
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		std::cout << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (pdb_end_time - pdb_start_time) << " seconds.\n";
		abinitio_protocol.clear_checkpoints();
	}
	jobdist.shutdown();

}

void ThisApplication::run() {
	setup();

	if (option [ rerun ]) {
		do_rerun();
		return;
	};
	fold();
	return;

}

void* my_main( void * )
{
	ThisApplication my_app;
	my_app.run();
	return NULL;
}

int
main( int argc, char * argv [] )
{
	ThisApplication::register_options( );
	//FoldConstraint::register_options()
	init( argc, argv );
	protocols::viewer::viewer_main( my_main );
	return 0;
}


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
#define DIM6 6
#define XX 0
#define YY 1
#define ZZ 2
static inline void oprod(const rvec a,const rvec b,rvec c)
{
	c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
	c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
	c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

void jacobi(double a[6][6],double d[],double v[6][6],int *nrot)
{
	int j,i;
	int iq,ip;
	double tresh,theta,tau,t,sm,s,h,g,c;
	double b[DIM6];
	double z[DIM6];
	int const n( DIM6 );
	for (ip=0; ip<n; ip++) {
		for (iq=0; iq<n; iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=0; ip<n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1; i<=50; i++) {
		sm=0.0;
		for (ip=0; ip<n-1; ip++) {
			for (iq=ip+1; iq<n; iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=0; ip<n-1; ip++) {
			for (iq=ip+1; iq<n; iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
						&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=0; j<ip; j++) {
						ROTATE(a,j,ip,j,iq)
		}
					for (j=ip+1; j<iq; j++) {
						ROTATE(a,ip,j,j,iq)
						}
					for (j=iq+1; j<n; j++) {
						ROTATE(a,ip,j,iq,j)
						}
					for (j=0; j<n; j++) {
						ROTATE(v,j,ip,j,iq)
						}
					++(*nrot);
				}
			}
		}
		for (ip=0; ip<n; ip++) {
			b[ip] +=  z[ip];
			d[ip]  =  b[ip];
			z[ip]  =  0.0;
		}
	}
	assert(0);
}



void calc_fit_R(int natoms,rvec *xp,rvec *x,matrix R)
{

	int    c,r,n,j,i,irot;
	double omega[ DIM6 ][ DIM6 ];
	double om[ DIM6 ] [ DIM6 ];
	double d[ DIM6 ],xnr,xpc;
	matrix vh,vk,u;
	Real   mn;
	int    index;
	Real   max_d;

	for(i=0; i<DIM6; i++) {
		d[i]=0;
		for(j=0; j<DIM6; j++) {
			omega[i][j]=0;
			om[i][j]=0;
		}
	}

	/* clear matrix U */
	for ( int i=0; i<DIM;i++)
		for ( int j=0; j<DIM; j++) u[i][j]=0;
	/*calculate the matrix U*/
	for(n=0;(n<natoms);n++) {
		if ((mn = 1.0) != 0.0) {
			for(c=0; (c<DIM); c++) {
				xpc=xp[n][c];
				for(r=0; (r<DIM); r++) {
					xnr=x[n][r];
					u[c][r]+=mn*xnr*xpc;
				}
			}
		}
	}
	dump_matrix(DIM, u);
	/*construct omega*/
	/*omega is symmetric -> omega==omega' */
	for(r=0; r<DIM6; r++)
		for(c=0; c<=r; c++)
			if (r>=DIM && c<DIM) {
				omega[r][c]=u[r-DIM][c];
				omega[c][r]=u[r-DIM][c];
			} else {
				omega[r][c]=0;
				omega[c][r]=0;
			}
	dump_matrix(DIM6, omega);
	/*determine h and k*/
	jacobi( omega,d,om,&irot);
	/*real   **omega = input matrix a[0..n-1][0..n-1] must be symmetric
	 *int     natoms = number of rows and columns
	 *real      NULL = d[0]..d[n-1] are the eigenvalues of a[][]
	 *real       **v = v[0..n-1][0..n-1] contains the vectors in columns
	 *int      *irot = number of jacobi rotations
	 */
	dump_matrix( 2*DIM, omega );
	dump_matrix ( 2*DIM, om );
	index=0; /* For the compiler only */

	/* Copy only the first two eigenvectors */
	for(j=0; j<2; j++) {
		max_d=-1000;
		for(i=0; i<DIM6; i++)
			if (d[i]>max_d) {
				max_d=d[i];
				index=i;
			}
		d[index]=-10000;
		for(i=0; i<DIM; i++) {
			vh[j][i]=sqrt(2.0)*om[i][index];
			vk[j][i]=sqrt(2.0)*om[i+DIM][index];
		}
	}
	/* Calculate the last eigenvector as the outer-product of the first two.
	 * This insures that the conformation is not mirrored and
	 * prevents problems with completely flat reference structures.
	 */

	dump_matrix( DIM, vh );
	dump_matrix( DIM, vk );
	oprod(vh[0],vh[1],vh[2]);
	oprod(vk[0],vk[1],vk[2]);
	dump_matrix( DIM, vh );
	dump_matrix( DIM, vk );

	/*determine R*/
	for(r=0; r<DIM; r++)
		for(c=0; c<DIM; c++)
			R[r][c] = vk[0][r]*vh[0][c] +
					vk[1][r]*vh[1][c] +
					vk[2][r]*vh[2][c];
	dump_matrix( DIM, R );
}


#if 0
	fill_coordinates( pose, ifit_, x );

	FArray1D< numeric::Real > ww( nfit_, 1.0 );
	FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;
	trDebug << "fit group \n";
	for ( int i = 1; i<=nfit_ ; i++) {
		for ( int d = 1; d<=3; d++ ) {
			trDebug << x(d, i) << " ";
		}
		trDebug << "\n";
	}
	trDebug << "--------------------\n";
	//find rotation translation structure
	numeric::rms::findUU( x, xref_, ww, nfit_, uu, ctx );

	for ( int i = 1; i<=nfit_ ; i++) {
		for ( int d = 1; d<=3; d++ ) {
			trDebug << x(d, i) << " ";
		}
		trDebug << "\n";
	}
	trDebug << "--------------------\n";

	//	fill_coordinates( pose, ipca_, x );

// 	// align center of mass to origin
// 	for ( int k = 1; k <= 3; ++k ) {
// 		Real temp1 = 0.0;
// 		for ( Size j = 1; j <= npca_; ++j ) {
// 			temp1 += x(k,j);
// 		}
// 		temp1 /= 1.0*npca_;
// 		for ( Size j = 1; j <= npca_; ++j ) {
// 			x(k,j) -= temp1;
// 		}
// 	}

// 	for ( int i=1; i<=npca_; i++) {
// 		FArray1D< numeric::Real > xdum( 3, 0.0 );
// 		for ( int d=1; d<=3; d++ ) {
// 			xdum ( d ) = 0;
// 			for ( int k=1; k<=3; k++ ) {
// 				xdum( d ) += uu( d, k) * x( k, i );
// 			}
// 		}
// 		for ( int d=1; d<=3; d++ ) {
// 			x( d, i) = xdum( d);// - xav_(d, i);
// 		}
// 	}
	trDebug << "rotated and translated:\n" ;
	for ( int i = 1; i<=npca_ ; i++) {
		for ( int d = 1; d<=3; d++ ) {
			trDebug << x(d, i) << " ";
		}
		trDebug << "\n";
	}

	//Compute projection
	proj.resize( nvec_ );
	for ( Size v = 1; v <= nvec_; v++ ) {
		proj[ v ]=0;
		for ( Size k = 1; k <= npca_; k++ ) {
			for ( Size d = 1; d <= 3; d++ ) {
				proj[ v ]+= x( d, k)*eigvec_( d, k, v)/10.0;
			}
		}
	}

#endif
