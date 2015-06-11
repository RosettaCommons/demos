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

import re
import warnings

from PDB.PDBIO import PDBIO


_hydrogen=re.compile("[123 ]*H.*")


class ChainSelector:
    """
    Only accepts residues with right chainid
    and between start and end. Remove hydrogens, waters and ligands.
    Only use model 0 by default.
    """
    def __init__(self, chain_id, start, end, model_id=0):
        self.chain_id=chain_id
        self.start=start
        self.end=end
        self.model_id=0

    def accept_model(self, model):
        # model - only keep model 0
        if model.get_id()==self.model_id:
            return 1
        return 0

    def accept_chain(self, chain):
        if chain.get_id()==self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        # residue - between start and end
        hetatm_flag, resseq, icode=residue.get_id()
        if hetatm_flag!=" ":
            # skip HETATMS
            return 0
        if icode!=" ":
            warnings.warn("WARNING: Icode at %s" % residue.get_id(),
                          RuntimeWarning)
        if self.start<=resseq<=self.end:
            return 1
        return 0

    def accept_atom(self, atom):
        # atoms - get rid of hydrogens
        name=atom.get_id()
        if _hydrogen.match(name):
            return 0
        else:
            return 1


def extract(structure, chain_id, start, end, filename):
    """
    Write out selected portion to filename.
    """
    sel=ChainSelector(chain_id, start, end)
    io=PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)



if __name__=="__main__":

    from PDB.PDBParser import PDBParser

    import sys

    p=PDBParser()
    s=p.get_structure("scr", sys.argv[1])

    extract(s, " ", 1, 100, "out.pdb")

