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
import shutil
from os import popen, system, mkdir
from os.path import basename, exists, isdir
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
import Util
from glob import glob

class PrepareRound2:
    def __init__( self, options ):
        self.opts_ = options

        self.r1_mc_picked_frags_dict_   = self.get_fragfiles_from_a_dir( self.opts_.round1_mc_picked_frags_dir )
        self.r1_picked_pos_             = self.r1_mc_picked_frags_dict_.keys()

        self.r1_frags_dict_             = self.get_fragfiles_from_a_dir( self.opts_.round1_frags_dir )
        self.r1_placement_results_dict_ = self.parse_placement_results_file( self.opts_.round1_clustering_selected_results_file )
        self.r1_placed_pos_             = self.r1_frags_dict_.keys()

        self.r2_frags_dict_             = self.get_fragfiles_from_a_dir( self.opts_.round2_frags_dir )
        self.r2_placement_results_dict_ = self.parse_placement_results_file( self.opts_.round2_clustering_selected_results_file )
        self.r2_placed_pos_             = self.r2_frags_dict_.keys()

        if options.total_rsd_list:
            self.total_rsd_positions_   = [ int(l.strip()) for l in open( options.total_rsd_list, "r" ).readlines() ]
        elif options.total_rsd:
            self.total_rsd_positions_   = range( 1, options.total_rsd+1 )
        else:
            self.total_rsd_positions_   = range( 1, max( self.r1_placement_results_dict_.keys() )+1 )
        print "total residue positions: %s" %( " ".join( map( str, self.total_rsd_positions_ ) ))
        self.round1_postions_needed_for_round2()



    def fragfn_to_pos( self, fragfile ):
        pos = Util.extract_fragid( fragfile )[1]
        return pos # it's int


    def get_fragfiles_from_a_dir( self, dir ):
        ''' store frag filename from a dir
            dict[ position ] = [ fragfile1, fragfile2, ... ] '''

        dict = {}
        assert isdir( dir ), dir

        fragfiles = map( basename, glob( dir+"/after*pdb" ) )

        assert len( fragfiles ) >0

        for fragfile in fragfiles:
            #print fragfile
            pos = self.fragfn_to_pos( fragfile )

            if dict.has_key( pos ):
                dict[ pos ].append( fragfile )
            else:
                dict[ pos ] = [ fragfile ]

        return dict


    def parse_placement_results_file( self, results_file ):
        ''' store placement results file as
            dict[ position ] = lines '''

        dict = {}
        with open( results_file, "r" ) as f:
            for l in f:
                if l.startswith("#"): continue
                ls = l.strip().split()
                pos = int(ls[0])
                if dict.has_key( pos ):
                    dict[ pos ] += l
                else:
                    dict[ pos ] = l
        return dict


    def round1_postions_needed_for_round2( self ):
        ''' some positions from r1 have no fragments being assigned
            even though they have some other frags covered that region
            > total_pos - r1_assigned_pos - r2_placed_pos '''

        # make sure they don't have any residues in common
        common_rsds = list( set(self.r1_picked_pos_) & set(self.r2_placed_pos_) )
        assert ( len( common_rsds ) == 0 ), "rsd: " + " ".join( map( str, common_rsds ) ) + " in common"

        self.round1_positions_needed_for_round2_ = list( set(self.total_rsd_positions_) - set(self.r1_picked_pos_) - set(self.r2_placed_pos_) )


    def prepare_fragfiles( self ):

        ''' 1. it copoes frags '''

        if not exists( self.opts_.dest_dir ):
            mkdir( self.opts_.dest_dir )
        else:
            system("rm -rf %s" % self.opts_.dest_dir )
            mkdir( self.opts_.dest_dir )

        # copy unassigned frags from round1
        for pos in self.round1_positions_needed_for_round2_:
            try:
                fragfiles = self.r1_frags_dict_[ pos ]
            except:
                stderr.write("WARNING: position %s has no placements from round1. Something is wrong.\n" % pos )
                continue

            for fragfile in self.r1_frags_dict_[ pos ]:
                fragfile = self.opts_.round1_frags_dir + "/" + fragfile
                #print fragfile
                assert exists( fragfile )
                shutil.copy( fragfile, self.opts_.dest_dir )

        # copy assigned frags from round1
        system("cp %s/* %s" %( self.opts_.round1_mc_picked_frags_dir, self.opts_.dest_dir ))

        # copy frags from round2
        system("cp %s/* %s" %( self.opts_.round2_frags_dir, self.opts_.dest_dir ))


    def dump_all_clustering_selected_file( self ):
        '''parsed position that havent assigned from r1_mc_picked and r2_run from round1_clustering_selected_file
        '''
        outlines = ""
        for pos in self.total_rsd_positions_:

            if pos in self.round1_positions_needed_for_round2_:
                try:
                    outlines += self.r1_placement_results_dict_[ pos ]
                    continue
                except:
                    stderr.write("WARNING: position %s has no resultlines from round1. Something is wrong.\n" % pos )
                    continue

            elif pos in self.r2_placed_pos_:
                try:
                    outlines += self.r2_placement_results_dict_[ pos ]
                    continue
                except:
                    stderr.write("WARNING: position %s has no resultlines from round2. Something is wrong.\n" % pos )
                    continue

            elif pos in self.r1_picked_pos_:
                picked_fragfile = self.r1_mc_picked_frags_dict_[ pos ][0]
                for l in self.r1_placement_results_dict_[ pos ].split("\n"):
                    if not l.strip(): continue
                    if l.startswith("#"): continue
                    ls = l.strip().split()
                    fragfile = ls[4]
                    if picked_fragfile == fragfile:
                        outlines += l+("\n")
                        continue

            else:
                stderr.write("ERROR: something is really wrong here\n")
                exit()

        out = open( self.opts_.outfile, "w" )
        out.write( outlines )
        out.close()




if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-a",  "--round1_mc_picked_frags_dir", required=True, help="")
    parser.add_argument("-f1", "--round1_clustering_selected_results_file", required=True, help="")
    parser.add_argument("-r1", "--round1_frags_dir", required=True, help="")

    parser.add_argument("-f2", "--round2_clustering_selected_results_file", required=True, help="")
    parser.add_argument("-r2", "--round2_frags_dir", required=True, help="")

    parser.add_argument("-d",  "--dest_dir", default="combined_round1_round2_frags", help="output dir")
    parser.add_argument("-o",  "--outfile", default="combined_round1_round2_frags.sc", help="")
    parser.add_argument("-l",  "--total_rsd_list", help="")
    parser.add_argument("-t",  "--total_rsd", type=int, help="")
    args = parser.parse_args()


    round2 = PrepareRound2( args )
    round2.prepare_fragfiles()
    round2.dump_all_clustering_selected_file()
