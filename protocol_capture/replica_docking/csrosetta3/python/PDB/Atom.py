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

"""Atom class, used in Structure objects."""
try:
   import numpy
except ImportError:
   pass 

from PDB.Entity import DisorderedEntityWrapper
from PDB.PDBExceptions import PDBConstructionWarning
from PDB.Vector import Vector
from PDB import IUPACData

class Atom:
    def __init__(self, name, coord, bfactor, occupancy, altloc, fullname, serial_number,
                 element=None):
        """
        Atom object.

        The Atom object stores atom name (both with and without spaces),
        coordinates, B factor, occupancy, alternative location specifier
        and (optionally) anisotropic B factor and standard deviations of
        B factor and positions.

        @param name: atom name (eg. "CA"). Note that spaces are normally stripped.
        @type name: string

        @param coord: atomic coordinates (x,y,z)
        @type coord: Numeric array (Float0, size 3)

        @param bfactor: isotropic B factor
        @type bfactor: number

        @param occupancy: occupancy (0.0-1.0)
        @type occupancy: number

        @param altloc: alternative location specifier for disordered atoms
        @type altloc: string

        @param fullname: full atom name, including spaces, e.g. " CA ". Normally
        these spaces are stripped from the atom name.
        @type fullname: string

        @param element: atom element, e.g. "C" for Carbon, "HG" for mercury,
        @type fullname: uppercase string (or None if unknown)

        """
        self.level="A"
        # Reference to the residue
        self.parent=None
        # the atomic data
        self.name=name      # eg. CA, spaces are removed from atom name
        self.fullname=fullname  # e.g. " CA ", spaces included
        self.coord=coord
        self.bfactor=bfactor
        self.occupancy=occupancy
        self.altloc=altloc
        self.full_id=None   # (structure id, model id, chain id, residue id, atom id)
        self.id=name        # id of atom is the atom name (e.g. "CA")
        self.disordered_flag=0
        self.anisou_array=None
        self.siguij_array=None
        self.sigatm_array=None
        self.serial_number=serial_number
        # Dictionary that keeps addictional properties
        self.xtra={}
        assert not element or element == element.upper(), element
        self.element = self._assign_element(element)
        self.mass = self._assign_atom_mass()

    def _assign_element(self, element):
        """Tries to guess element from atom name if not recognised."""
        if not element or element.capitalize() not in IUPACData.atom_weights:
            # Inorganic elements have their name shifted left by one position
            #  (is a convention in PDB, but not part of the standard).
            # isdigit() check on last two characters to avoid mis-assignment of
            # hydrogens atoms (GLN HE21 for example)

            if self.fullname[0] != " " and not self.fullname[2:].isdigit():
                putative_element = self.name.strip()
            else:
                # Hs may have digit in [0]
                if self.name[0].isdigit():
                    putative_element = self.name[1]
                else:
                    putative_element = self.name[0]

            if putative_element.capitalize() in IUPACData.atom_weights:
                msg = "Used element %r for Atom (name=%s) with given element %r" \
                      % (putative_element, self.name, element)
                element = putative_element
            else:
                msg = "Could not assign element %r for Atom (name=%s) with given element %r" \
                      % (putative_element, self.name, element)
                element = ""
            import warnings
            warnings.warn(msg, PDBConstructionWarning)

        return element

    def _assign_atom_mass(self):
        # Needed for Bio/Struct/Geometry.py C.O.M. function
        if self.element:
            return IUPACData.atom_weights[self.element.capitalize()]
        else:
            return float('NaN')


    # Special methods

    def __repr__(self):
        "Print Atom object as <Atom atom_name>."
        return "<Atom %s>" % self.get_id()

    def __sub__(self, other):
        """
        Calculate distance between two atoms.

        Example:
            >>> distance=atom1-atom2

        @param other: the other atom
        @type other: L{Atom}
        """
        diff=self.coord-other.coord
        return numpy.sqrt(numpy.dot(diff,diff))

    # set methods

    def set_serial_number(self, n):
        self.serial_number=n

    def set_bfactor(self, bfactor):
        self.bfactor=bfactor

    def set_coord(self, coord):
        self.coord=coord

    def set_altloc(self, altloc):
        self.altloc=altloc

    def set_occupancy(self, occupancy):
        self.occupancy=occupancy

    def set_sigatm(self, sigatm_array):
        """
        Set standard deviation of atomic parameters.

        The standard deviation of atomic parameters consists
        of 3 positional, 1 B factor and 1 occupancy standard
        deviation.

        @param sigatm_array: standard deviations of atomic parameters.
        @type sigatm_array: Numeric array (length 5)
        """
        self.sigatm_array=sigatm_array

    def set_siguij(self, siguij_array):
        """
        Set standard deviations of anisotropic temperature factors.

        @param siguij_array: standard deviations of anisotropic temperature factors.
        @type siguij_array: Numeric array (length 6)
        """
        self.siguij_array=siguij_array

    def set_anisou(self, anisou_array):
        """
        Set anisotropic B factor.

        @param anisou_array: anisotropic B factor.
        @type anisou_array: Numeric array (length 6)
        """
        self.anisou_array=anisou_array


    # Public methods

    def flag_disorder(self):
        """Set the disordered flag to 1.

        The disordered flag indicates whether the atom is disordered or not.
        """
        self.disordered_flag=1

    def is_disordered(self):
        "Return the disordered flag (1 if disordered, 0 otherwise)."
        return self.disordered_flag

    def set_parent(self, parent):
        """Set the parent residue.

        Arguments:
        o parent - Residue object
        """
        self.parent=parent

    def detach_parent(self):
        "Remove reference to parent."
        self.parent=None

    def get_sigatm(self):
        "Return standard deviation of atomic parameters."
        return self.sigatm_array

    def get_siguij(self):
        "Return standard deviations of anisotropic temperature factors."
        return self.siguij_array

    def get_anisou(self):
        "Return anisotropic B factor."
        return self.anisou_array

    def get_parent(self):
        "Return parent residue."
        return self.parent

    def get_serial_number(self):
        return self.serial_number

    def get_name(self):
        "Return atom name."
        return self.name

    def get_id(self):
        "Return the id of the atom (which is its atom name)."
        return self.id

    def get_full_id(self):
        """Return the full id of the atom.

        The full id of an atom is the tuple
        (structure id, model id, chain id, residue id, atom name, altloc).
        """
        return self.parent.get_full_id()+((self.name, self.altloc),)

    def get_coord(self):
        "Return atomic coordinates."
        return self.coord

    def get_bfactor(self):
        "Return B factor."
        return self.bfactor

    def get_occupancy(self):
        "Return occupancy."
        return self.occupancy

    def get_fullname(self):
        "Return the atom name, including leading and trailing spaces."
        return self.fullname

    def get_altloc(self):
        "Return alternative location specifier."
        return self.altloc

    def get_level(self):
        return self.level

    def transform(self, rot, tran):
        """
        Apply rotation and translation to the atomic coordinates.

        Example:
                >>> rotation=rotmat(pi, Vector(1,0,0))
                >>> translation=array((0,0,1), 'f')
                >>> atom.transform(rotation, translation)

        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        self.coord=numpy.dot(self.coord, rot)+tran

    def get_vector(self):
        """
        Return coordinates as Vector.

        @return: coordinates as 3D vector
        @rtype: Vector
        """
        x,y,z=self.coord
        return Vector(x,y,z)


class DisorderedAtom(DisorderedEntityWrapper):
    """
    This class contains all Atom objects that represent the same disordered
    atom. One of these atoms is "selected" and all method calls not caught
    by DisorderedAtom are forwarded to the selected Atom object. In that way, a
    DisorderedAtom behaves exactly like a normal Atom. By default, the selected
    Atom object represents the Atom object with the highest occupancy, but a
    different Atom object can be selected by using the disordered_select(altloc)
    method.
    """
    def __init__(self, id):
        """
        Arguments:
        o id - string, atom name
        """
        self.last_occupancy=-1
        DisorderedEntityWrapper.__init__(self, id)

    # Special methods

    def __repr__(self):
        return "<Disordered Atom %s>" % self.get_id()

    def disordered_add(self, atom):
        "Add a disordered atom."
        # Add atom to dict, use altloc as key
        atom.flag_disorder()
        # set the residue parent of the added atom
        residue=self.get_parent()
        atom.set_parent(residue)
        altloc=atom.get_altloc()
        occupancy=atom.get_occupancy()
        self[altloc]=atom
        if occupancy>self.last_occupancy:
            self.last_occupancy=occupancy
            self.disordered_select(altloc)
