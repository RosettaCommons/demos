# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from rosetta import *
import math
import sys

class DOF(object):
  minval = None
  maxval = None
  def __init__(self,pose,seqpos):
    self.pose = pose
    self.seqpos = seqpos

  def set_dof_value(self,val):
#    print >>sys.stderr, type(self), type(val),type(min),type(max)
    if self.maxval: val = float(min(val,self.maxval))
    if self.minval: val = float(max(val,self.minval))
    self._set_pose_dof(val)

  def __iadd__(self,incr):
    val = self._get_pose_dof()
    self.set_dof_value( val+incr )

  def _set_pose_dof(self,val):
    raise NotImplementedError,"set_pose_dof is abstract in class DOF!"

  def _get_pose_dof(self):
    raise NotImplementedError,"set_pose_dof is abstract in class DOF!"


class PhiDOF(DOF):
  minval = 0.0
  maxval = 360.0
  def _get_pose_dof(self):
#    print >>sys.stderr, type(self),'get pose dof ',self.pose,self.seqpos
    return self.pose.phi(self.seqpos)

  def _set_pose_dof(self,val):
#    print >>sys.stderr, type(self), type(val), 'setting pose dof ',self.pose,self.seqpos,val
    self.pose.set_phi(self.seqpos,val)



