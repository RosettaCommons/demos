# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *
sys.stdout = sys.stderr

init("/Users/sheffler/svn/branches/minirosetta_database")

pose  = Pose(pdb_file="input/1rdg.pdb")

score = ScoreFunction( default=True )

class PyMover(Mover):
  count = 0
  def apply(self,pose,move_map):
    self.count += 1
    print "PyMover called %i times," % self.count,
    print "score is %f\n!" % score(pose)

py_mover = PyMover( name="pymover", weight=10.0 )

relax = ClassicRelax( score, pose )

relax.moves1().add_move( new_move=py_mover )

relax.run()



