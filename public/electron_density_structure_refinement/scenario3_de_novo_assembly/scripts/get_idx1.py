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
from ScoreTable import ScoreTable
from sys import stderr, exit, stdout
import cPickle as pickle
import os


if __name__=='__main__':
    '''
    converting indexed scorefiles into npy
    '''
    parser = ArgumentParser()
    parser.add_argument("-p", "--selected_frags_path", required=True, type=str, help="")
    parser.add_argument("-c", "--clustering_results", required=True, type=str, help="")
    args = parser.parse_args()


    # fragidx file
    fragidx_pickle = args.selected_frags_path + "/index_Dict.pickle"
    assert os.path.exists( fragidx_pickle )
    frag_idx1 = open( "frags.idx1", "w" )
    dict = pickle.load( open( fragidx_pickle, "rb" ) )[1]
    for idx in dict.keys():
        frag_idx1.write("%s %s %s %s %s\n" %( idx, dict[idx][0], dict[idx][1], dict[idx][2], dict[idx][3] ))
    frag_idx1.close()


    # get density_score_Dict.pickle"
    scoretable = ScoreTable( args.selected_frags_path )
    scoretable.density_score_reader( args.clustering_results )

    # output all_density_rmsd.idx1
    den_idx1 = open( "all_density.idx1", "w" )
    dict = pickle.load( open( "density_score_Dict.pickle", "rb" ) )
    for pos in dict.keys():
        for idx in dict[pos]:
            if idx < 0: continue
            den_idx1.write("%s %s %s\n" %( idx, dict[pos][idx][0], dict[pos][idx][1] ))



