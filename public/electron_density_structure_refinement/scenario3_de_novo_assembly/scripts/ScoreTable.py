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
import numpy as np
from glob import glob
from pprint import pprint
from os.path import basename, exists, abspath
from Util import *
from sys import stderr, stdout
import cPickle as pickle
from Containers import FragID


class ScoreTable:
    '''
    *** caveat: the frag index using for the cleaner has to be the same
    This class encapsulates:
        0. indexing dict:        self._frag_index_dict
        1. one-body score table: self._density_score_dict
        2. two-body score table: self._score_table
    '''

    def __init__( self, selected_frags_path ):
        self._frag_to_index_dict = {}
        self._index_to_frag_dict = {}

        self.indexing( selected_frags_path ) # get indexing dictionaries
        print self.__dimension


    def retrieve_density_score_dict( self ):
        return self._density_score_dict


    def indexing( self, selected_frags_path  ):
        '''
        _frag_to_index_dict[ frag_id_tuple ] = index
        _index_to_frag_dict[ index ] = frag_id_tuple
        _score_table: initialize a 3D array:
        '''

        pickle_fn = "index_Dict.pickle"
        selected_frags_path = abspath( selected_frags_path ) + "/"

        if exists( selected_frags_path + pickle_fn ):
            stdout.write("Loading %s ..." % pickle_fn ); stdout.flush()
            pkl = open( selected_frags_path + pickle_fn, "rb" )
            ( self._frag_to_index_dict, self._index_to_frag_dict ) = pickle.load( pkl ); stdout.write(" Done!\n")

        else:
            stdout.write("Loading frag_fns ..."); stdout.flush()
            frag_list = map( basename, glob( selected_frags_path + "after*pdb" ))

            sorted_frags = sorted( frag_list ) # sorted in order to make the indexing determinisitic

            # is going to change the index to be 1-somenumber
            # 140221 for vector1 starting with 1
            #for index in range( 1, len( sorted_frags )+1 ):
            for index in range( len( sorted_frags ) ):
                frag_fn = sorted_frags[ index ]
                frag_id_tuple = extract_fragid( frag_fn )

                self._frag_to_index_dict[ frag_id_tuple ] = index+1
                self._index_to_frag_dict[ index+1 ] = frag_id_tuple

            pickle.dump( ( self._frag_to_index_dict, self._index_to_frag_dict ), open( selected_frags_path + pickle_fn, "w") ); stdout.write(" Done!\n")

        self.__dimension = len( self._frag_to_index_dict ) # score table dimension ( square )

        #self.initialize_twobody_score_table()


    def initialize_twobody_score_table( self ):
        size = len( self._frag_to_index_dict )
        #self._twobody_score_table = np.zeros(( size, size, 3 ))  # 3 <- 0: overlap_score; 1: ( closab_score, clash_score )
        self._overlapping_check_table = np.zeros(( size, size ))  # 3 <- 0: overlap_score; 1: ( closab_score, clash_score )

        self._overlap_score_table     = np.zeros(( size, size ))  # 3 <- 0: overlap_score; 1: ( closab_score, clash_score )
        self._clash_score_table       = np.zeros(( size, size ))  # 3 <- 0: overlap_score; 1: ( closab_score, clash_score )
        self._closab_score_table      = np.zeros(( size, size ))  # 3 <- 0: overlap_score; 1: ( closab_score, clash_score )


    def get_candidate_frags( self, pos ):
        ''' This is going to be called at MonteCarlo object through ScoreFunction abject.
            Return frag_idxs for a given residue position
        '''
        assert isinstance( pos, int )
        assert self._density_score_dict.has_key( pos )
        return self._density_score_dict[ pos ]


    def get_density_score_dict( self ):
        ''' return a dictionary that:
            { pos: { frag_idx : ( densityScore, rmsd ) } }
            see also - get_tuple_density_score_dict()
        '''
        return self._density_score_dict


    def get_tuple_density_score_dict( self ):
        ''' The density_score_dict that:
            { pos: { fragid_tuple : ( densityScore, rmsd ) } }
            see also - get_density_score_dict()
        '''
        density_score_dict = {}

        for pos in self._density_score_dict.keys():
            for each_frag_idx in self._density_score_dict[ pos ]:
                data = self._density_score_dict[ pos ][ each_frag_idx ]
                if pos in density_score_dict.keys():
                    density_score_dict[ pos ][ self.index_to_frag( each_frag_idx ) ] = data
                else:
                    density_score_dict[ pos ] = { self.index_to_frag( each_frag_idx ) : data }

        return density_score_dict


    def idx_to_fragpos( self, frag_idx ):
        ''' given a frag_idx, return the position of it '''
        return self.index_to_frag( frag_idx )[1]


    def index_to_frag( self, index ):
        ''' return a fragid tuple '''
        if index < 0:
            return ( 9, -index, 0, 0 )
        else:
            return self._index_to_frag_dict[ index ]


    def frag_to_index( self, frag ):
        if isNullFrag( FragID( frag ) ):
            return -frag[1]
        else:
            return self._frag_to_index_dict[ frag ]


    def rmsd_lookup( self, frag_idx ):
        return self._density_score_dict[ self.idx_to_fragpos( frag_idx ) ][ frag_idx ][1]


    def density_score_lookup( self, frag_idx ):
        return self._density_score_dict[ self.idx_to_fragpos( frag_idx ) ][ frag_idx ][0]


    def overlapping_check_table_lookup( self, ifrag_idx, jfrag_idx ):
        ''' precalculate overlapping_check '''
        return self._overlapping_check_table[ ifrag_idx, jfrag_idx ]


    def overlap_score_lookup( self, ifrag_idx, jfrag_idx ):
        return self._overlap_score_table[ ifrag_idx, jfrag_idx ]


    def closab_score_lookup( self, ifrag_idx, jfrag_idx ):
        return self._closab_score_table[ ifrag_idx, jfrag_idx ]


    def clash_score_lookup( self, ifrag_idx, jfrag_idx ):
        return self._clash_score_table[ ifrag_idx, jfrag_idx ]


    def score_files_reader( self, density_scorefile, overlap_scorefile, nonoverlap_scorefile ):
        '''
        '''
        self.density_score_reader( density_scorefile )
        self.overlap_score_reader( overlap_scorefile )
        self.nonoverlap_score_reader( nonoverlap_scorefile )


    def overlapping_check( self, i_fragidx, j_fragidx ):
        ''' > This is arguably to be one of the most important function in the class
            1. Check whether two fragments are overlapping or not
            2. If one fragments is overlapping, only the overlap score (overlap/clash) is being calculated, not nonoverlap score (closab/clash)
        '''
        i_frag_id = self.index_to_frag( i_fragidx )
        j_frag_id = self.index_to_frag( j_fragidx )

        i_mer, i_pos, j_mer, j_pos = i_frag_id[0], i_frag_id[1], j_frag_id[0], j_frag_id[1]

        # swap_poses, want to make sure j_pos is always further than i_pos
        if ( i_pos > j_pos ):
            j_mer, j_pos, i_mer, i_pos = i_mer, i_pos, j_mer, j_pos

        assert( j_pos > i_pos ), "i_pos: %s; j_pos: %s" %( i_pos, j_pos )
        offset = j_pos - i_pos  #should therefore always be positive

        # "=" means there is one-residue overlapping
        if offset <= i_mer - 1:
            return True
        else:
            return False


    def make_square_array( self, size ):
        '''
        just want to reduce typing ( would it? )
        '''
        return np.zeros(( size, size ))


    def overlap_score_reader( self, scorefile ):
        ''' '''
        npy_fn = "overlap_score.npy"
        if exists( npy_fn ):
            stdout.write("Loading %s... " % npy_fn ); stdout.flush()
            self._overlap_score_table = np.load( npy_fn ); stdout.write("Done!\n")
        else:
            self._overlap_score_table = self.make_square_array( self.__dimension )

            with open( scorefile, "r" ) as f:
                stdout.write("Reading %s...\n" % scorefile ); stdout.flush()
                for l in f:
                    if l.startswith("#"): continue
                    ls = l.split()
                    assert len(ls)==3
                    ifrag_idx      = int( ls[0] )
                    jfrag_idx      = int( ls[1] )
                    overlap_score = float( ls[2] )

                    self._overlap_score_table[ ifrag_idx, jfrag_idx ] = overlap_score
                    self._overlap_score_table[ jfrag_idx, ifrag_idx ] = overlap_score

            self._overlap_score_table.dump("overlap_score.npy")


    def nonoverlap_score_reader( self, scorefile ):
        ''' '''
        if exists("closab_score.npy") and exists("clash_score.npy") and exists("overlapping_check_table.npy"):

            stdout.write("Loading overlapping_check_table.npy... "); stdout.flush()
            self._overlapping_check_table = np.load("overlapping_check_table.npy"); stdout.write("Done!\n")

            stdout.write("Loading closab_score.npy... " ); stdout.flush()
            self._closab_score_table = np.load("closab_score.npy"); stdout.write("Done!\n")

            stdout.write("Loading clash_score.npy... " ); stdout.flush()
            self._clash_score_table = np.load("clash_score.npy"); stdout.write("Done!\n")

        else:
            self._overlapping_check_table = self.make_square_array( self.__dimension )
            self._closab_score_table      = self.make_square_array( self.__dimension )
            self._clash_score_table       = self.make_square_array( self.__dimension )

            with open( scorefile, "r" ) as f:
                stdout.write("Reading %s...\n" % scorefile ); stdout.flush()
                for l in f:
                    if l.startswith("#"): continue
                    ls = l.split()
                    assert len(ls)==4
                    ifrag_idx     = int( ls[0] )
                    jfrag_idx     = int( ls[1] )
                    closab_score  = float( ls[2] )
                    clash_score   = int( ls[3] )
                    overlapping   = self.overlapping_check( ifrag_idx, jfrag_idx )

                    self._overlapping_check_table[ ifrag_idx, jfrag_idx ] = overlapping
                    self._overlapping_check_table[ jfrag_idx, ifrag_idx ] = overlapping

                    self._closab_score_table[ ifrag_idx, jfrag_idx ] = closab_score
                    self._closab_score_table[ jfrag_idx, ifrag_idx ] = closab_score

                    self._clash_score_table[ ifrag_idx, jfrag_idx ] = clash_score
                    self._clash_score_table[ jfrag_idx, ifrag_idx ] = clash_score


            self._overlapping_check_table.dump("overlapping_check_table.npy")
            self._closab_score_table.dump("closab_score.npy")
            self._clash_score_table.dump("clash_score.npy")


    def density_score_reader( self, scorefile ):
        ''' '''
        pickle_fn = "density_score_Dict.pickle"
        self._density_score_dict = {}

        if exists( pickle_fn ):
            stdout.write("Loading %s... " % pickle_fn )
            pkl = open( pickle_fn, "rb" )
            self._density_score_dict = pickle.load( pkl ); stdout.write("Done!\n")

        else:
            stdout.write( "Reading %s ..." % scorefile ); stdout.flush()
            self._density_score_dict = {} # ( rank, mer ) : densityScore

            with open( scorefile, "r" ) as f:
                for l in f:
                    ls = l.split()
                    assert len( ls ) <= 6, l # tag for round2 round1

                    densityScore = float( ls[1] )
                    rmsd         = float( ls[2] )

                    fragid_tuple = extract_fragid( ls[4] )
                    frag_idx = self._frag_to_index_dict[ fragid_tuple ]
                    pos = fragid_tuple[1]

                    if pos in self._density_score_dict.keys():
                        self._density_score_dict[ pos ][ frag_idx ] = ( densityScore, rmsd )
                    else:
                        # for the first-time initialization
                        self._density_score_dict[ pos ] = { frag_idx : ( densityScore, rmsd ) }
                        self._density_score_dict[ pos ][ -pos ] = ( 0, 0 ) # for each position, create a Null fragment fake density score, because you will need to print out density score when you select a null frag for a position

            pickle.dump( self._density_score_dict, open("density_score_Dict.pickle", "w") ); stdout.write("done\n")


    ## to prepare idx1 idx2 scores
    def raw_overlap_score_cleaner( self, overlap_scorefile ):
        assert exists( overlap_scorefile )

        stderr.write( "reading %s\n" % overlap_scorefile )
        outlines = ""
        with open( overlap_scorefile, "r" ) as f:
             for l in f:
                if l.startswith("#"): continue
                ls = l.split()
                assert len(ls)==5
                ifrag_idx     = self._frag_to_index_dict[ extract_fragid( ls[0] ) ]
                jfrag_idx     = self._frag_to_index_dict[ extract_fragid( ls[1] ) ]
                overlap_score = float( ls[3] )
                outlines += "%s %s %s\n" %( ifrag_idx, jfrag_idx, overlap_score )

        return outlines


    def raw_nonoverlap_score_cleaner( self, nonoverlap_scorefile, closab_weights_dict=None, nonclosab_weight=1 ):
        ''' for closab score, by calculating gap size you can decide db to use '''
        assert exists( nonoverlap_scorefile )
        stderr.write( "reading %s\n" % nonoverlap_scorefile )
        if closab_weights_dict:
            assert isinstance( closab_weights_dict, dict )

        outlines = ""
        with open( nonoverlap_scorefile, "r" ) as f:
            for l in f:
                if l.startswith("#"): continue
                if "fragA" in l: continue
                ls = l.split()
                assert len(ls)==8
                ifrag_idx     = self._frag_to_index_dict[ extract_fragid( ls[0] ) ]
                jfrag_idx     = self._frag_to_index_dict[ extract_fragid( ls[1] ) ]
                gap_size      = int( float( ls[3] ) )
                closab_score  = float( ls[5] )

                if closab_score == -1: # only reweight closab ones
                    if closab_weights_dict:
                        if closab_weights_dict.has_key( gap_size ): # gap_size fall in between 1 - 25
                            closab_score  = closab_score*closab_weights_dict[ gap_size ]
                        else:
                            stderr.write("ERROR: this shouldn't happen - weight dict doesn't have gap_size %s \n" % gap_size )

                elif closab_score == 1: # Nonclosable - almost 99% is correct for all the gap sizes
                    closab_score  = closab_score*nonclosab_weight

                clash_score   = int( float( ls[6] ) )
                outlines += "%s %s %s %s\n" %( ifrag_idx, jfrag_idx, closab_score, clash_score )

        return outlines


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument("-f", "--frags", nargs="+", required=True, help="")
    parser.add_argument("-o", "--overlap_scorefile", required=True, help="")
    parser.add_argument("-n", "--nonoverlap_scorefile", required=True, help="")
    args = parser.parse_args()

    #outfile = open( args.scorefile+".clean.idx", "w" )
    #lines = ""

    scoretable = ScoreTable( [ basename( frag ) for frag in args.frags ] )
    scoretable.score_file_reader( args.overlap_scorefile, args.nonoverlap_scorefile )

    #if args.scoretype == "overlap":
    #    lines = scoretable.raw_overlap_score_cleaner( args.scorefile )
    #elif args.scoretype == "nonoverlap":
    #    lines = scoretable.raw_nonoverlap_score_cleaner( args.scorefile )

    #outfile.write( lines )
    #outfile.close()




