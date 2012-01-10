# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *
init("/Users/sheffler/svn/branches/minirosetta_database")
pose = Pose()
pose_from_pdb(pose,"input/1rdg.pdb")
score = ScoreFunction()
score.set('fa_atr',1)
score(pose)
score.set('fa_atr',0.8)
score(pose)
score.set('fa_rep',0.8)
score(pose)
score.set('fa_sol',0.8)
score(pose)
score(pose)
score.set('fa_sol',0.7)
score(pose)
rlx = ClassicRelax()
rlx = ClassicRelax(score,pose)

for t,w in weights: s.set(t,w)
