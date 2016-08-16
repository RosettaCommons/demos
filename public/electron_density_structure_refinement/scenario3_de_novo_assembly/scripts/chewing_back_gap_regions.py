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
from os import popen, system
from os.path import basename, exists
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
import denovo_model_building_util

def get_chunks( seq_dict ):
    ''' given a seq_dict[ rsd ] = state (1: assigned; 0: unassigned)
        return a chunk_dict[ idx ] = [ continuous chunk residue numbers ]
    '''
    chunk_dict = {}
    print seq_dict
    i=1
    idx=""
    while i <= len(seq_dict.keys()):
        while seq_dict[i] == 1:
            if not idx:
                idx=i
            if not chunk_dict.has_key( idx ):
                chunk_dict[ idx ] = [i]
            else:
                chunk_dict[ idx ].append(i)
            i+=1
            print idx, chunk_dict[idx]
        i+=1
        idx=i
    return chunk_dict


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-p", "--pdb", required=True, help="")
    parser.add_argument("-f", "--fasta", required=True, help="")
    parser.add_argument("-c", "--chewing", default=4, type=int, help="")
    args = parser.parse_args()

    seq = ( denovo_model_building_util.fasta_file_reader( args.fasta ) )
    seq_dict = {}
    for i in range(1,len(seq)+1):
        seq_dict[i]=0

    xyz_dict, junk, pdbline_dict, k = denovo_model_building_util.create_xyzDict( args.pdb )
    for rsd in xyz_dict.keys():
        if seq_dict.has_key( rsd ):
            seq_dict[rsd] = 1


    chunk_dict = get_chunks(seq_dict)
    for chunk in sorted( chunk_dict.keys() ):
        #print len(chunk_dict[chunk]), chunk_dict[chunk]
        #print "REMARK", len(chunk_dict[chunk][args.chewing:-args.chewing]), chunk_dict[chunk][args.chewing:-args.chewing]
        for rsd in chunk_dict[chunk][args.chewing:-args.chewing]:
            stdout.write( pdbline_dict[ rsd ] )

