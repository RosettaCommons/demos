# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *
from numeric import *
from DOF import *
from random import *
from math import *

init(dbstr)

class PyMover(Mover):
  def __init__(self,name,dofs,magnitude=1,weight=1):
    self.magnitude = magnitude
    self.dofs = dofs
    Mover.__init__(self,name,weight)

  def apply(self,pose,move_map):
    print >>sys.stderr,"PyMover(mag=%f) ignoring MoveMap and using its own DOFs!"%self.magnitude
    for dof in self.dofs:
      # print >>sys.stderr,"PyMover DOF:",dof
      dof += (2*random()-1)*self.magnitude

class PyScore(ScoreFunction):
  def __call__(self,pose):
    print >>sys.stderr, "USING PyScore!!! looping over atoms...slow..."
    s = ScoreFunction.__call__(self,pose)
    # for at1,at2 in gen_heavy_atom_pairs(pose):
    #   if distance_squared( at1.xyz(), at2.xyz() ) < 16:
    #     s -= 0.02
    return s

pose  = Pose(pdb_file="input/1rdg.pdb")
score = PyScore( default=True )
print score(pose)

mm = MoveMap()
mm.set_bb(True)
opt = MinimizerOptions("linmin",1,True)
mzr = AtomTreeMinimizer()
mzr.run(pose,mm,score,opt)

dofs = [PhiDOF(pose,seqpos) for seqpos in range(1,pose.total_residue())]
py_mover = PyMover( name="pymover", dofs=dofs, magnitude=1,weight=1000000.0 )

relax = ClassicRelax( score, pose )
# relax.moves1().add_move( new_move=py_mover )
# relax.moves2().add_move( new_move=py_mover )
# relax.moves3().add_move( new_move=py_mover )
relax.run()




# STATUS reject      -112.086  -113.237  -114.804 shear
# STATUS low-accept  -114.996  -114.996  -116.460 shear





