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

"""Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Author: Thomas Hamelryck.  Additional code by Kristian Rother.
"""

# Get a Structure object from a PDB file
from PDBParser import PDBParser

try:
    # Get a Structure object from an mmCIF file
    from MMCIFParser import MMCIFParser
except:
    # Not compiled I guess 
    pass

# Download from the PDB
from PDBList import PDBList 

# Parse PDB header directly
from parse_pdb_header import parse_pdb_header

# Find connected polypeptides in a Structure
from Polypeptide import PPBuilder, CaPPBuilder, is_aa, standard_aa_names
# This is also useful :-)
from PDB.Raf import to_one_letter_code

# IO of PDB files (including flexible selective output)
from PDBIO import PDBIO, Select

# Some methods to eg. get a list of Residues
# from a list of Atoms.
import Selection

# Superimpose atom sets
from Superimposer import Superimposer

# 3D vector class
from Vector import Vector, calc_angle, calc_dihedral, refmat, rotmat, rotaxis,\
        vector_to_axis, m2rotaxis, rotaxis2m

# Alignment module
from StructureAlignment import StructureAlignment

# DSSP handle 
# (secondary structure and solvent accessible area calculation)
from DSSP import DSSP, make_dssp_dict

# Residue depth: 
# distance of residue atoms from solvent accessible surface
from ResidueDepth import ResidueDepth, get_surface

# Calculation of Half Sphere Solvent Exposure
from HSExposure import HSExposureCA, HSExposureCB, ExposureCN

# Kolodny et al.'s backbone libraries
from FragmentMapper import FragmentMapper

# Write out chain(start-end) to PDB file
from Dice import extract

# Fast atom neighbor search
# Depends on KDTree C++ module
try:
    from NeighborSearch import NeighborSearch
except ImportError:
    pass
