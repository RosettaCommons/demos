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

"""Map the residues of two structures to each other based on a FASTA alignment
file.
"""

from PDB.Raf import to_one_letter_code

from PDB import Selection
from PDB.Polypeptide import is_aa


class StructureAlignment:
    """
    This class aligns two structures based on an alignment of their
    sequences.
    """
    def __init__(self, fasta_align, m1, m2, si=0, sj=1):
        """
        fasta_align --- Alignment object 
        m1, m2 --- two models
        si, sj --- the sequences in the Alignment object that
                correspond to the structures
        """
        l=fasta_align.get_alignment_length()
        s1=fasta_align.get_seq_by_num(si)
        s2=fasta_align.get_seq_by_num(sj)
        # Get the residues in the models
        rl1=Selection.unfold_entities(m1, 'R')
        rl2=Selection.unfold_entities(m2, 'R')
        # Residue positions
        p1=0
        p2=0
        # Map equivalent residues to each other
        map12={}
        map21={}
        # List of residue pairs (None if -)
        duos=[]
        for i in range(0, l):
            column=fasta_align.get_column(i)
            aa1=column[si]
            aa2=column[sj]
            if aa1!="-":
                # Position in seq1 is not -
                while 1:
                    # Loop until an aa is found
                    r1=rl1[p1]
                    p1=p1+1
                    if is_aa(r1):
                        break
                self._test_equivalence(r1, aa1)
            else:
                r1=None
            if aa2!="-":
                # Position in seq2 is not -
                while 1:
                    # Loop until an aa is found
                    r2=rl2[p2]
                    p2=p2+1
                    if is_aa(r2):
                        break
                self._test_equivalence(r2, aa2)
            else:
                r2=None
            if r1:
                # Map residue in seq1 to its equivalent in seq2
                map12[r1]=r2
            if r2:
                # Map residue in seq2 to its equivalent in seq1
                map21[r2]=r1
            # Append aligned pair (r is None if gap)
            duos.append((r1, r2))
        self.map12=map12
        self.map21=map21
        self.duos=duos

    def _test_equivalence(self, r1, aa1):
        "Test if aa in sequence fits aa in structure."
        resname=r1.get_resname()
        resname=to_one_letter_code[resname]
        assert(aa1==resname)

    def get_maps(self):
        """
        Return two dictionaries that map a residue in one structure to 
        the equivealent residue in the other structure.
        """
        return self.map12, self.map21

    def get_iterator(self):
        """
        Iterator over all residue pairs.
        """
        for i in range(0, len(self.duos)):
            yield self.duos[i]


if __name__=="__main__":
    import sys
    from PDB.Alphabet import generic_protein
    #from Bio import AlignIO
    from PDB import PDBParser

    if len(sys.argv) != 4:
        print "Expects three arguments,"
        print " - FASTA alignment filename (expect two sequences)"
        print " - PDB file one"
        print " - PDB file two"
        sys.exit()

    # The alignment
    fa=AlignIO.read(open(sys.argv[1]), "fasta", generic_protein)

    pdb_file1=sys.argv[2]
    pdb_file2=sys.argv[3]

    # The structures
    p=PDBParser()
    s1=p.get_structure('1', pdb_file1)
    p=PDBParser()
    s2=p.get_structure('2', pdb_file2)

    # Get the models
    m1=s1[0]
    m2=s2[0]

    al=StructureAlignment(fa, m1, m2)

    # Print aligned pairs (r is None if gap)
    for (r1,r2) in al.get_iterator():
        print r1, r2

