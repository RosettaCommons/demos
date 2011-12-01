# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
from _numeric import *

def print_matrix(self):
  return "%(xx)f %(xy)f %(xz)f\n%(yx)f %(yy)f %(yz)f\n%(zx)f %(zy)f %(zz)f\n"%self

def print_vector(self):
  return "%(x)f %(y)f %(z)f"%self

strfunc = { "xyzMatrix":print_matrix, "xyzVector":print_vector }

def getitem(self,k):
  """docstring for getitem"""
  return getattr(self,k)

for n in dir():
  types = ('xyzVector','xyzMatrix')
  for t in types:
    if n.startswith("%s"%t):
      getattr(__import__(__name__),n).__getitem__ = getitem
      getattr(__import__(__name__),n).__str__  = strfunc[t]
      getattr(__import__(__name__),n).__repr__ = strfunc[t]

xyzVector_double.__mul__ = lambda u,v: xyzVector_double( u.x*v.x , u.y*v.y , u.z*v.z )


def test():
  if 'xyzVector_double' in dir():
    v = xyzVector_double(1,2,3)
    x_unit = xyzVector_double(1,0,0)
    y_unit = xyzVector_double(0,1,0)
    z_unit = xyzVector_double(0,0,1)
    equal_length(v,v)
    print min(v, 0.5*v)
    print max(3*v,v)
    M = xyzMatrix_double()
    R = rotation_matrix( x_unit, 3.1415/2 ) # rotate around x axis
    P = projection_matrix( y_unit ) #project along y axis
    T = P*R # rotate then project
    print T*v

if __name__ == '__main__':
  test()

