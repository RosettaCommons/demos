#! /usr/bin/python
# Copyright (c) 2010 Robert L. Campbell
#
# A PyMOL script for drawing a CGO plane from the coordinates of three atoms (pk1,pk2,pk3 by default)

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#
# @brief Script to visualize membrane planes in pymol
# @author JKLeman (julia.koehler1982@gmail.com)

from pymol.cgo import *
from pymol import cmd

# Store xyz coordinates
class XYZCoord():
	
  def __init__(self, x=0, y=0, z=0):
    self.x = x
    self.y = y
    self.z = z
		
  def __str__(self):
    return "(%8.3f, %8.3f, %8.3f)" % (self.x, self.y, self.z)

cmd.extend( "XYZCoord", XYZCoord )

################################################################################

# Subtract v2 from v1, represented as XYZcoord classes (see above)
def subtract( v1, v2 ):

  x = v1.x - v2.x
  y = v1.y - v2.y
  z = v1.z - v2.z
				
  vf = XYZCoord( x, y, z)
  return vf

cmd.extend( "subtract", subtract )

################################################################################

# Add v1 and v2
def add( v1, v2 ):

  x = v1.x + v2.x
  y = v1.y + v2.y
  z = v1.z + v2.z
				
  vf = XYZCoord( x, y, z)
  return vf

cmd.extend( "add", add )

################################################################################

# Cross-product
def cross( v1, v2 ):

  x = v1.y * v2.z - v1.z * v2.y
  y = v1.z * v2.x - v1.x * v2.z
  z = v1.x * v2.y - v1.y * v2.x
				
  vf = XYZCoord( x, y, z)
  return vf

cmd.extend( "cross", cross )

################################################################################

# Return length of vector v, represented as XYZcoord class (see above)
def length( v ):

  return math.sqrt( v.x**2 + v.y**2 + v.z**2 )

cmd.extend( "length", length )

################################################################################

# Normalize vector v to length l
def normalize( v, l ):

  v_length = length( v )

  x = v.x * l / v_length
  y = v.y * l / v_length
  z = v.z * l / v_length

  vf = XYZCoord( x, y, z)

  return vf

cmd.extend( "normalize", normalize )

################################################################################

# A helper function for computing the normal to a triangular facet
def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):

  nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)
  ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)
  nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)

  return (nx,ny,nz)



def draw_plane_cgo(name,apex1,apex2,apex3,apex4,color,transparency):

  """
DESCRIPTION
    Create a CGO plane from three arbitary coordinates

USAGE
    draw_plane_cgo apex1, apex2, apex3, apex4, color

    where each apex is a 3-element vector and color is a 3-element RGB
    list defining the color of the plane (where each value of R, G
    and B is between 0 and 1 inclusive).

  """

  # Convert args to floating point numbers
  x1,y1,z1 = map(float,apex1)
  x2,y2,z2 = map(float,apex2)
  x3,y3,z3 = map(float,apex3)
  x4,y4,z4 = map(float,apex4)
  if type(color) == type(''):
    color = map(float,color.replace('(','').replace(')','').split(','))

  # Compute the normal vector for the triangle
  normal1 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)
  normal2 = compute_normal(x1, y1, z1, x3, y3, z3, x4, y4, z4)
  normal3 = compute_normal(x2, y2, z2, x3, y3, z3, x4, y4, z4)

  # Create the CGO objects
  obj = [

    BEGIN, TRIANGLE_STRIP,
	
#    COLOR, color[0], color[1], color[2],
	ALPHA, 0.5,
	COLOR, 0.5, 0.5, 0.5,
    NORMAL, normal1[0], normal1[1], normal1[2],
    VERTEX, x1, y1, z1,
    VERTEX, x2, y2, z2,
    VERTEX, x3, y3, z3,
    VERTEX, x4, y4, z4,

    END
  ]

  # Display them
  cmd.load_cgo(obj,name)

def draw_plane(name,atom1='(pk1)',atom2='(pk2)',atom3='(pk3)',atom4='(pk4)',color=[0,0,0],transparency=0.5):
  """
DESCRIPTION
    Create a CGO plane from four atomic coordinates

USAGE
    draw_plane name, atom1, atom2, atom3, atom4, color

    where each atom is a standard PyMOL selection (defaults to pk1,pk2
    and pk3) and color is a 3-element RGB tuple defining the color
    of the plane (where each value of R, G and B is between 0 and 1
    inclusive).  The color defaults to (1,1,1).

  """

# get coordinates from atom selections
  coor1 = cmd.get_model(atom1).atom[0].coord
  coor2 = cmd.get_model(atom2).atom[0].coord
  coor3 = cmd.get_model(atom3).atom[0].coord
  coor4 = cmd.get_model(atom4).atom[0].coord
  draw_plane_cgo(name,coor1,coor2,coor3,coor4,color,transparency)

cmd.extend("draw_plane", draw_plane)
