#!/usr/bin/env python
import os,shutil
from optparse import OptionParser
from set_alignment import *
from write_pdb_files import *

sa = set_alignment()
wpf = write_pdb_files()

"""
Makes an alignment to a zinc site in a protein. Atom names of
ligand has to be set in main in the list SUBSTRATEATOMS

@param -f: pdbfile with zincprotein
@type f  : string
@param -l: pdbfile with ligand
@type l  : string
@return  : aligned_coordinates in file named 'aligned_ligand.pdb'

"""
# The atom names should be orderly given 
# with zinc ion first and then het atom
# It is assumed that Metal is located before
# Heteroatom in pdb file.

SUBSTRATEATOMS = ['ZN1','O2']

# Constants for output file
METAL = 'ZN'


def main():
    parser = OptionParser()
    parser.add_option('-f',dest='pdbfile',
                      help='Protein file with zinc ion')
    parser.add_option('-l',dest='ligfile',
                      help='PDB format file with ligand')
    (options,args) = parser.parse_args()
    # String of atom names from the ligand to be aligned with
    # Set the path to the pdb file
    PDB = options.pdbfile
    LIGANDFILE = options.ligfile
    # other methods into set_alignment
    sa.generate_alignment_all(METAL,PDB,LIGANDFILE,SUBSTRATEATOMS)
    # Add hack for correcting for the right ZN xyz
    zn_line = wpf.get_metal_ion(PDB)
    
    c_file = wpf.set_correct_line('sample_substrate.pdb',zn_line[0])
    wpf.write_aligned_coordinates(c_file,'aligned_ligand.pdb')

    # Clean files
    tmp_files =['cry.pdb','tmp.pdb','sample_substrate.pdb','superimposed.pdb']
    for i in tmp_files:
        os.remove(i)
    
if __name__ == "__main__":
    main()
