# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *

class PyScore(ScoreFunction):
  def __call__(self,pose):
    s = ScoreFunction.__call__(self,pose)
    for atom1,atom2 in gen_atom_pairs(pose):
      if 4.0 > atom1.xyz().distance( atom2.xyz() ):
        s += 0.01
    return s

pose = Pose("input/1rdg.pdb")
my_score = PyScore(default=True)
my_relax = ClassicRelax( score_func=my_score, pose=pose )
my_relax.run()

sys.exit()




# (the utility function that
# generates Pose's atom pairs)
def gen_atom_pairs(pose):
  for res1 in pose:
    for res2 in pose:
      for atom1 in res1:
        for atom2 in res2:
          yield atom1,atom2

# this is how you loop over
# atom pairs in python:

for atom1,atom2 in gen_atom_pairs(pose):
  do_some_stuff( atom1, atom2 )

# this is how you loop over
# atom pairs in rosetta++:

sys.exit()



m_Temperature = 0.8;
lj_ramp_cycles=1;
lj_ramp_inner_cycles = 1;
start_rep_weight = 0.02;
end_rep_weight = 1.0;
stage2_repack_period = 1;
stage2_cycles = 1;
stage3_cycles = 1;

mm = MoveMap();
mm.set_bb ( True );
mm.set_chi( True );

small_mover = SmallMover()
small_mover.set_angle_max( 'H', 2.0 );
small_mover.set_angle_max( 'E', 2.0 );
small_mover.set_angle_max( 'L', 3.0 );
small_mover.nmoves( 5 );

shear_mover = ShearMover()
shear_mover.set_angle_max( 'H', 2.0 );
shear_mover.set_angle_max( 'E', 2.0 );
shear_mover.set_angle_max( 'L', 3.0 );
shear_mover.nmoves( 5 );



#print 'create packer task'
#task = create_packer_task( pose )

mc = MonteCarlo(pose,score,0.8)
mopt = MinimizerOptions("linmin",1.0,True)
#ropt = RottrialOptions( ptask, 0.01 )
#popt = ProtocolOptions( mopt, ropts)

print score( pose )
score.accumulate_residue_total_energies( pose )



