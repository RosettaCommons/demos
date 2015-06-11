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
# Copyright (C) 2006, Thomas Hamelryck (thamelry@binf.ku.dk)
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

"""Wrappers for PSEA, a program for secondary structure assignment.

See this citation for P-SEA, PMID: 9183534

Labesse G, Colloc'h N, Pothier J, Mornon J-P:  P-SEA: a new efficient
assignment of secondary structure from C_alpha.
Comput Appl Biosci 1997 , 13:291-295

ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/
"""

import os

from PDB.Polypeptide import is_aa


def run_psea(fname):
    """Run PSEA and return output filename.
    
    Note that this assumes the P-SEA binary is called "psea" and that it is
    on the path.
    
    Note that P-SEA will write an output file in the current directory using
    the input filename with extension ".sea".
    
    Note that P-SEA will write output to the terminal while run.
    """
    os.system("psea "+fname)
    last=fname.split("/")[-1]
    base=last.split(".")[0]
    return base+".sea"

def psea(pname):
    """Parse PSEA output file."""
    fname=run_psea(pname)
    start=0
    ss=""
    fp=open(fname, 'r')
    for l in fp.readlines():
        if l[0:6]==">p-sea":
            start=1
            continue
        if not start:
            continue
        if l[0]=="\n":
            break
        ss=ss+l[0:-1]
    fp.close()
    return ss

def psea2HEC(pseq):
    """Translate PSEA secondary structure string into HEC."""
    seq=[]
    for ss in pseq:
        if ss=="a":
            n="H"
        elif ss=="b":
            n="E"
        elif ss=="c":
            n="C"
        seq.append(n)
    return seq

def annotate(m, ss_seq):
    """Apply seconardary structure information to residues in model."""
    c=m.get_list()[0]
    all=c.get_list()
    residues=[]
    # Now remove HOH etc.
    for res in all:
        if is_aa(res):
            residues.append(res)
    L=len(residues)
    if not (L==len(ss_seq)):
        raise ValueError("Length mismatch %i %i" % (L, len(ss_seq)))
    for i in range(0, L):
        residues[i].xtra["SS_PSEA"]=ss_seq[i]
    #os.system("rm "+fname)

class PSEA:
    def __init__(self, model, filename):
        ss_seq=psea(filename)
        ss_seq=psea2HEC(ss_seq)
        annotate(model, ss_seq)
        self.ss_seq=ss_seq

    def get_seq(self):
        """
        Return secondary structure string.
        """
        return self.ss_seq
        

if __name__=="__main__":

    import sys
    from PDB import PDBParser

    # Parse PDB file
    p=PDBParser()
    s=p.get_structure('X', sys.argv[1])

    # Annotate structure with PSEA sceondary structure info
    PSEA(s[0], sys.argv[1])
