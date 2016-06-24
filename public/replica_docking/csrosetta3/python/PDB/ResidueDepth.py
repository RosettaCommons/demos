###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###
# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

#                 Biopython License Agreement
#
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
#
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE.

"""Calculation of residue depth using command line tool MSMS.

This module uses Michel Sanner's MSMS program for the surface calculation
(specifically commands msms and pdb_to_xyzr). See:
http://mgltools.scripps.edu/packages/MSMS

Residue depth is the average distance of the atoms of a residue from 
the solvent accessible surface.

Residue Depth:

    >>> rd = ResidueDepth(model, pdb_file)
    >>> print rd[(chain_id, res_id)]

Direct MSMS interface:

    Typical use:

        >>> surface = get_surface("1FAT.pdb")

    Surface is a Numeric array with all the surface 
    vertices.  

    Distance to surface:

        >>> dist = min_dist(coord, surface)

    where coord is the coord of an atom within the volume
    bound by the surface (ie. atom depth).

    To calculate the residue depth (average atom depth
    of the atoms in a residue):

        >>> rd = residue_depth(residue, surface)
"""

import os
import tempfile

try:
 import numpy
except:
 pass

from PDB import Selection
from PDB.AbstractPropertyMap import AbstractPropertyMap
from PDB.Polypeptide import is_aa


def _read_vertex_array(filename):
    """
    Read the vertex list into a Numeric array.
    """
    fp=open(filename, "r")
    vertex_list=[]
    for l in fp.readlines():
        sl=l.split()
        if not len(sl)==9:
            # skip header
            continue
        vl=map(float, sl[0:3])
        vertex_list.append(vl)
    fp.close()
    return numpy.array(vertex_list)

def get_surface(pdb_file, PDB_TO_XYZR="pdb_to_xyzr", MSMS="msms"):
    """
    Return a Numeric array that represents 
    the vertex list of the molecular surface.

    PDB_TO_XYZR --- pdb_to_xyzr executable (arg. to os.system)
    MSMS --- msms executable (arg. to os.system)
    """
    # extract xyz and set radii
    xyz_tmp=tempfile.mktemp()
    PDB_TO_XYZR=PDB_TO_XYZR+" %s > %s"
    make_xyz=PDB_TO_XYZR % (pdb_file, xyz_tmp)
    os.system(make_xyz)
    assert os.path.isfile(xyz_tmp), \
        "Failed to generate XYZR file using command:\n%s" % make_xyz
    # make surface
    surface_tmp=tempfile.mktemp()
    MSMS=MSMS+" -probe_radius 1.5 -if %s -of %s > "+tempfile.mktemp()
    make_surface=MSMS % (xyz_tmp, surface_tmp)
    os.system(make_surface)
    surface_file=surface_tmp+".vert"
    assert os.path.isfile(surface_file), \
        "Failed to generate surface file using command:\n%s" % make_surface
    # read surface vertices from vertex file
    surface=_read_vertex_array(surface_file)
    # clean up tmp files
    # ...this is dangerous
    #os.system("rm "+xyz_tmp)
    #os.system("rm "+surface_tmp+".vert")
    #os.system("rm "+surface_tmp+".face")
    return surface


def min_dist(coord, surface):
    """
    Return minimum distance between coord
    and surface.
    """
    d=surface-coord
    d2=numpy.sum(d*d, 1)
    return numpy.sqrt(min(d2))

def residue_depth(residue, surface):
    """
    Return average distance to surface for all
    atoms in a residue, ie. the residue depth.
    """
    atom_list=residue.get_unpacked_list()
    length=len(atom_list)
    d=0
    for atom in atom_list:
        coord=atom.get_coord()
        d=d+min_dist(coord, surface)
    return d/length

def ca_depth(residue, surface):
    if not residue.has_id("CA"):
        return None
    ca=residue["CA"]
    coord=ca.get_coord()
    return min_dist(coord, surface)

class ResidueDepth(AbstractPropertyMap):
    """
    Calculate residue and CA depth for all residues.
    """
    def __init__(self, model, pdb_file):
        depth_dict={}
        depth_list=[]
        depth_keys=[]
        # get_residue
        residue_list=Selection.unfold_entities(model, 'R')
        # make surface from PDB file
        surface=get_surface(pdb_file)
        # calculate rdepth for each residue
        for residue in residue_list:
            if not is_aa(residue):
                continue
            rd=residue_depth(residue, surface)
            ca_rd=ca_depth(residue, surface)
            # Get the key
            res_id=residue.get_id()
            chain_id=residue.get_parent().get_id()
            depth_dict[(chain_id, res_id)]=(rd, ca_rd)
            depth_list.append((residue, (rd, ca_rd)))
            depth_keys.append((chain_id, res_id))
            # Update xtra information
            residue.xtra['EXP_RD']=rd
            residue.xtra['EXP_RD_CA']=ca_rd
        AbstractPropertyMap.__init__(self, depth_dict, depth_keys, depth_list)


if __name__=="__main__":

    import sys
    from PDB import PDBParser

    p=PDBParser()
    s=p.get_structure("X", sys.argv[1])
    model=s[0]

    rd=ResidueDepth(model, sys.argv[1])


    for item in rd:
        print item

