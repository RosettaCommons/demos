#!/usr/bin/env python
import sys
from numpy import *

"""
Generate resfile for design - two inputs needed
@start > 0
@end


"""

def get_ligand_coordinates(pdbfile):
    pdbfl = open(pdbfile,'r')
    ligand_coor = []
    for line in pdbfl:
        if line[0:4] == 'HETA':
            x = str(line[31:39]).rstrip()
            y = str(line[39:47]).rstrip()
            z = str(line[47:55]).rstrip()
            ligand_coor.append(array([float(x),float(y),float(z)]))
    return ligand_coor

def get_length_protein(pdbfile):
    pdbfl = open(pdbfile,'r')
    first = '0'
    start = ''
    end = ''
    for line in pdbfl:
        if line[0:4] == 'ATOM':
            if first == '0':
                start = str(line[23:26]).rstrip()
                first = '1'
            elif first == '1':
                end = str(line[23:26]).rstrip()
    print start,end
    return int(start),int(end)
            

def get_residues_in_pdbfile(pdbfile):
    pdbfl = open(pdbfile,'r')
    tmp_residues = []
    for line in pdbfl:
        if line[0:4] == 'ATOM':
            tmp_residues.append(int(str(line[23:26]).rstrip()))
    unique_residues = []
    residues = set(tmp_residues)
    for i in residues:
        unique_residues.append(i)
    unique_residues.sort()
    return unique_residues


def main():
    # File name of pdb file
    pdbfile = sys.argv[1]
    # Collecting the ligand coordinates in
    # a vector
    ligand = get_ligand_coordinates(pdbfile)
    # Collecting residue id for residues within or outside
    # ligand
    within_ligand = []
    outside_ligand = []
    pdbfl = open(pdbfile,'r')
    for i in ligand:
        for line in pdbfl:
            if line[0:4] == 'ATOM':
                x = str(line[31:39]).rstrip()
                y = str(line[39:47]).rstrip()
                z = str(line[47:55]).rstrip()
                tmp_vector = array([float(x),float(y),float(z)])

                tmp_length = linalg.norm(tmp_vector - i)
                if tmp_length < 15:
                    pdb_resid = int(str(line[23:26]).rstrip())
                    within_ligand.append(pdb_resid)
    set(within_ligand)
    start,end = get_length_protein(pdbfile)
    # resfile
    rs_file = open('resfile','w')
    rs_file.write('start\n')


    unique_residues = get_residues_in_pdbfile(pdbfile)
    
    rs_1 = ' A NATAA\n'
    rs_2 = ' A NATRO\n'


    while start <= end:
#        print type(start)
        if start in within_ligand:
            rs_file.write('\t'+str(start)+rs_1)
        elif start in unique_residues:
            rs_file.write('\t'+str(start)+rs_2)
        start = start + 1

if __name__ == "__main__":
    main()
