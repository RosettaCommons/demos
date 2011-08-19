# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *

init("/Users/sheffler/svn/branches/minirosetta_database")

pose = Pose("input/1rdg.pdb")

score = ScoreFunction( fa_atr=0.80, fa_rep=0.44, n_ci_2b_score_types=0.65, fa_pair=0.49,
                       hbond_sr_bb=1.17, hbond_lr_bb=1.17, hbond_bb_sc=1.17, hbond_sc=1.10,
                       rama=0.20, fa_dun=0.56 )

class PyMover(Mover):
  count = 0
  def apply(self,pose,move_map):
    self.count += 1
    print "PyMover called %i times! score is now %f, dofs:"%(self.count,score(pose))
    for dof in move_map: print `dof`

m = PyMover("test")

my_relax = ClassicRelax(score,pose)
my_relax.moves1().add_move( m )
my_relax.run()





  for ( int i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
    conformation::Residue const & rsd( pose.residue( i ) );
    if ( rsd.is_protein() && mm.get_bb( i ) ) {
      char const ss( pose.secstruct( i ) );
      if ( angle_max.count( ss ) ) {
        Real const mx( angle_max.find( ss )->second );
        if ( mx > 0.0 ) {
          pos_list.push_back( std::make_pair( i, mx ) );
        }
      }
    }
  }

  // sanity
  if ( pos_list.empty() ) {
    std::cout << "no movable positions in small_moves!" << std::endl;
    return;
  }

  // get the rama potential
  scoring::Ramachandran const & rama
    ( scoring::ScoringManager::get_instance()->get_Ramachandran() );

  // how many moves to make
  int const num = std::max( Size(1), std::min( num_in, pos_list.size()/2 ) );

  // positions at which we've moved
  utility::vector1< int > already_moved;

  // now loop
  for ( int k=1; k<= num; ++k ) {
    int tries(0);
    while ( true ) {
      ++tries;
      if ( tries > 10000 ) {
        break;
      }
      std::pair< int, Real > const & p
        ( pos_list[ static_cast< int >( sm_RG.uniform() * pos_list.size() + 1 ) ] );
      int const j( p.first );
      Real const big_angle ( p.second );
      Real const small_angle ( big_angle/2.0 );

//  next three lines skip ends of structure !!
//  fix a logic error here: the closer to the end, the higher the probability
//  to skip it; and when it is 1 or total_residue, the probability should be
//  zero; and the probability distribution should be symmetrical for two ends
//  residue:    N-   1   2   3   4   5   6
//  residue:    C-   t  t-1 t-2 t-3 t-4 t-5
//  prob to skip:    1  0.8 0.6 0.4 0.2 0.0
//  -- chu

      // need to add this back, prob want a pose.distance_to_chain_end(i)
      // function
      //
      //end = total_residue - 5;
      //if ( j <= 5 && static_cast< int >(ran3()*5+1) >= j ) goto L401;
      //if ( j > end && static_cast< int >(ran3()*5) + end <= j ) goto L401;

      if ( std::find( already_moved.begin(), already_moved.end(), j ) !=
           already_moved.end() ) continue;

      conformation::Residue const & rsd( pose.residue(j) );

      Real const current_phi( pose.phi(j) );
      Real const current_psi( pose.psi(j) );

      Real const phi_tmp
        ( util::periodic_range( current_phi - small_angle + sm_RG.uniform() * big_angle,
                                360.0 ) );
      Real const psi_tmp
        ( util::periodic_range( current_psi - small_angle + sm_RG.uniform() * big_angle,
                                360.0 ) );

      Real const old_rama_score
        ( rama.eval_rama_score_residue( rsd.aa(), current_phi, current_psi ));

      Real const new_rama_score
        ( rama.eval_rama_score_residue( rsd.aa(), phi_tmp, psi_tmp ) );

      if ( new_rama_score > old_rama_score ) {
        Real const boltz_factor = (old_rama_score-new_rama_score)/temp;
        Real const probability = std::exp(std::max(Real(-40.0),boltz_factor) );
        if ( sm_RG.uniform() >= probability ) continue;
      }

      pose.set_phi( j, phi_tmp );
      pose.set_psi( j, psi_tmp );

      already_moved.push_back( j );
      break;
    } // while ( true )
  } // k=1,num
