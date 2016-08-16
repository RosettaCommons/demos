#!/usr/bin/env python2.7
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
from argparse import ArgumentParser
from multiprocessing import Pool
#from string import Template
from os.path import exists, basename
from os import system, getcwd, popen
from sys import stderr, stdout, exit
import shutil

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

def pdb2fasta(pdb_fn):
    lines = open(pdb_fn,'r').readlines()

    fastaid = ""

    oldresnum = '   '
    count = 0;
    for line in lines:
        if (len(line)>20): # and (chainid == line[21]):
            line_edit = line
            if line[0:3] == 'TER':
                fastaid += '\n'
                #    break
            elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit)>75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]

            if line_edit[0:4] == 'ATOM':
                resnum = line_edit[23:26]
                if not resnum == oldresnum:
                    count = count + 1
                    longname = line_edit[17:20]
                    if longer_names.has_key(longname):
                        fastaid += longer_names[longname]
                    else:
                        fastaid += 'X'
                oldresnum = resnum

                newnum = '%3d' % count
                line_edit = line_edit[0:23] + newnum + line_edit[26:]

    return fastaid.strip()


def fasta_file_reader( fasta_fn ):
    from Bio import SeqIO
    fasta_seq = str( SeqIO.read( fasta_fn, "fasta" ).seq )
    return fasta_seq


def offset_pdb( pdb, startnum ):
    assert exists( pdb )
    out_pdblines = ""
    checked_startnum_in_pdb = False
    with open( pdb, "r" ) as f:
        for line in f:
            if line.startswith( "ATOM" ):
                rsn      = int( line[22:26]) # res number
                if not checked_startnum_in_pdb:
                    assert rsn <= startnum, "%s %s" %(rsn, startnum)
                    offset = startnum - rsn
                    checked_startnum_in_pdb = True
                    if offset == 0:
                        stderr.write("WARNING: offset = 0\n")

                rsn += offset
                out_pdblines += line[:22] + "%4s" % rsn + line[26:]

    return out_pdblines


def reindex_Nterminal_truncated_pdb( truncated_pdb_fn, fasta_fn ):
    pdb_seq = pdb2fasta( truncated_pdb_fn )
    fasta_seq = fasta_file_reader( fasta_fn )
    idx = fasta_seq.index( pdb_seq )

    stdout.write("detect N-terminus trimmed, shift the numbering to start at %s\n" %( idx + 1 ))
    stdout.write(fasta_seq + "\n")
    stdout.write("-"*idx + pdb_seq + "\n" )
    start_num = idx+1
    return start_num




if __name__=='__main__':
    parser = ArgumentParser( description="Important: This script can only be used to deals with the ''''terminus'''' truncated models - it will renumber the index truncated pdb" )
    parser.add_argument("-f", "--full_length_fasta_fn", required=True, help="")
    parser.add_argument("-i", "--missing_rsd_pdbs", required=True, nargs="+", help="")
    options = parser.parse_args()

    cwd = getcwd()
    for pdb in options.missing_rsd_pdbs:
        assert exists( pdb )
        startnum = reindex_Nterminal_truncated_pdb( pdb, options.full_length_fasta_fn )
        pdblines = offset_pdb( pdb, startnum )
        pdb_fn = basename(pdb).split(".pdb")[0] + ".renum.pdb"
        outstream = open( pdb_fn, "w" )
        outstream.write( pdblines )
        outstream.close()



