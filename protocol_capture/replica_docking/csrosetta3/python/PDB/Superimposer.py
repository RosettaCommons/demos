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

"""Superimpose two structures."""
try:
  import numpy
except ImportError:
  pass

#from Bio.SVDSuperimposer import SVDSuperimposer
#from PDB.PDBExceptions import PDBException


class Superimposer:
    """
    Rotate/translate one set of atoms on top of another,
    thereby minimizing the RMSD.
    """
    def __init__(self):
        self.rotran=None
        self.rms=None

    def set_atoms(self, fixed, moving):
        """
        Put (translate/rotate) the atoms in fixed on the atoms in 
        moving, in such a way that the RMSD is minimized.

        @param fixed: list of (fixed) atoms
        @param moving: list of (moving) atoms 
        @type fixed,moving: [L{Atom}, L{Atom},...]
        """
        if not (len(fixed)==len(moving)):
            raise PDBException("Fixed and moving atom lists differ in size")
        l=len(fixed)
        fixed_coord=numpy.zeros((l, 3))
        moving_coord=numpy.zeros((l, 3))
        for i in range(0, len(fixed)):
            fixed_coord[i]=fixed[i].get_coord()
            moving_coord[i]=moving[i].get_coord()
        sup=SVDSuperimposer()
        sup.set(fixed_coord, moving_coord)
        sup.run()
        self.rms=sup.get_rms()
        self.rotran=sup.get_rotran()

    def apply(self, atom_list):
        """
        Rotate/translate a list of atoms.
        """
        if self.rotran is None:
            raise PDBException("No transformation has been calculated yet")
        rot, tran=self.rotran
        rot=rot.astype('f')
        tran=tran.astype('f')
        for atom in atom_list:
            atom.transform(rot, tran)


if __name__=="__main__":
    import sys

    from PDB import PDBParser, Selection

    p=PDBParser()
    s1=p.get_structure("FIXED", sys.argv[1])
    fixed=Selection.unfold_entities(s1, "A")

    s2=p.get_structure("MOVING", sys.argv[1])
    moving=Selection.unfold_entities(s2, "A")

    rot=numpy.identity(3).astype('f')
    tran=numpy.array((1.0, 2.0, 3.0), 'f')

    for atom in moving:
        atom.transform(rot, tran)
    
    sup=Superimposer()

    sup.set_atoms(fixed, moving)

    print sup.rotran
    print sup.rms

    sup.apply(moving)
