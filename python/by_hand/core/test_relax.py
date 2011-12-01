# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
import sys
sys.path.append('lib')
from _test_relax import *


# class Shake(Mover):
#   def apply(self,pose,move_map):
#     for dof in move_map:
#       print dof

init("/Users/sheffler/svn/branches/minirosetta_database")

mm = MoveMap()
mm.set_bb(10,True)
print mm.get_bb(10)
print mm.get_bb(5)



#p = Pose()
#pose_from_pdb(p,"test_in.pdb")
#s = ScoreFunction()

#for t,w in [("fa_atr",0.80),("fa_rep",0.44),("fa_sol",0.65),("fa_pair",0.49),("hbond_bb_sc",1.17), \
#            ("fa_dun",0.56),("rama",0.2),("hbond_lr_bb",1.17),("hbond_sr_bb",1.17),("hbond_sc",1.10)]:
#  s.set_weight(getattr(ScoreType,t),w)

#my_relax = ClassicRelax(s,p)
#my_relax.run()

