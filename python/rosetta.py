# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
import sys
sys.path.append('by_hand/core/lib')
from _test_relax import *
from _test import *

dbstr = "/Users/sheffler/svn/branches/minirosetta_database"

DO_NOT_ERASE = []


def replace_init(cls,init):
  cls.__old_init = cls.__init__
  cls.__init__ = init


def py_clone(self):
  pe = type(self)()
  DO_NOT_ERASE.append(pe)
  return pe

################# ScoreFunction ######################

def set_weight(score,score_type,weight,pose=None):
  if not score_type in ScoreType.values.values():
    if not hasattr(ScoreType,score_type):
      print "WARNING weight not set: no score type %s!!!"%score_type
      if True: return None
    score_type = getattr( ScoreType, score_type )
  if not score_type in ScoreType.values.values():
    print "WARNING weight not set: no score type %s!!!"%score_type
    if True: return None
  oldscore = newscore = None
  if pose: oldscore = score(pose)
  score.set_weight( score_type , weight )
  if pose: newscore = score(pose)
  if pose: print 'old score: %f, new score: %f, delta: '%(oldscore,newscore,newscore-oldscore)
  return weight

def get_weight(score,score_type):
  if not hasattr(ScoreType,score_type):
    #print "WARNING weight doesn't exist: no score type %s!!!"%score_type
    if True: return
  t = getattr( ScoreType, score_type )
  return score.get_weight(t)

ScoreFunction.set_weight_py = set_weight
ScoreFunction.get_weight_py = get_weight
ScoreFunction.__getattr__ = get_weight
ScoreFunction.__setattr__ = set_weight

def scorefunc_to_str(score,joiner=', '):
  s = []
  for t in ScoreType.values.values():
    w = score.get_weight(t)
    if w > 0: s.append("%s=%.2f"%(t.name,w))
  return 'ScoreFunction( ' + joiner.join(s)+" )"

ScoreFunction.__repr__ = scorefunc_to_str
ScoreFunction.__str__  = lambda s: scorefunc_to_str(s,"\n               ")

def new_score_init(score,default=False,**args):
  ScoreFunction.__old_init(score)
  if default:
    score.fa_atr  = 0.80; score.hbond_sr_bb = 1.17;
    score.fa_pair = 0.49; score.hbond_bb_sc = 1.17;
    score.fa_rep  = 0.44; score.hbond_lr_bb = 1.17;
    score.fa_dun  = 0.56; score.hbond_sc    = 1.10;
    score.rama    = 0.20; score.n_ci_2b_score_types = 0.65;
  for t,w in args.items(): score.set_weight_py(t,w)

replace_init(ScoreFunction,new_score_init)

###################################
def new_relax_init(self,score_func,pose):
  self.__old_init(score_func,pose)

replace_init(ClassicRelax,new_relax_init)


#################### Atom #######################
def atom_to_str(a):
  return "Atom( type="+str(a.type()) +", coords:"+str(a.xyz())

Atom.__str__  = atom_to_str
Atom.__repr__ = atom_to_str


################ Resdiue ######################
def res_to_str(r):
  return "Residue( aa=%s, natom=%i )"%(r.name3(),r.natoms())

Residue.__str__  = res_to_str
Residue.__repr__ = res_to_str

##################### Pose ##################################
def new_pose_init(pose,pdb_file=None):
  Pose.__old_init(pose)
  pose.input_filename = pdb_file
  # pose._dofs = []
  if pdb_file: pose_from_pdb(pose,pdb_file)

replace_init(Pose,new_pose_init)
# Pose.dofs = lambda pose: pose._dofs

def pose_to_str(pose):
  return "Pose('%s'), nres=%i"%(pose.input_filename,pose.total_residue())

Pose.__repr__ = pose_to_str

##########################################
def new_mover_init(self,name,weight=1):
  Mover.__old_init(self,name)
  self.set_weight(weight)

replace_init(Mover,new_mover_init)

##############################################

ContextIndependentOneBodyEnergy.indicate_required_context_graphs = lambda x,y: None
ContextIndependentOneBodyEnergy.atomic_interaction_cutoff = lambda x: 0.0
ContextIndependentOneBodyEnergy.clone = py_clone


def gen_atom_pairs(pose,dist=None):
  for res1 in pose:
    for res2 in pose:
      for atom1 in res1:
        for atom2 in res2:
          if dist:
            if dist > numeric.distance( atom1.xyz(), atom2.xyz() ): continue
          yield atom1,atom2

def gen_heavy_atom_pairs(pose,dist=None):
  for res1 in pose:
    for atom1 in res1:
      if atom1.type() > 21: continue
      for res2 in pose:
        for atom2 in res2:
          if atom1.type() > 21: continue
          if dist:
            if dist > numeric.distance( atom1.xyz(), atom2.xyz() ): continue
          yield atom1,atom2


