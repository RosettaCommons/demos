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
from os.path import basename, exists
from os import getcwd

# from distance count
gap_weights_dict = { 1: 1.0000,
                     2: 0.5563,
                     3: 0.3529,
                     4: 0.2907,
                     5: 0.1735,
                     6: 0.1149,
                     7: 0.1032,
                     8: 0.0580,
                     9: 0.0469,
                    10: 0.0545,
                    11: 0.0457,
                    12: 0.0504,
                    13: 0.0510,
                    14: 0.0416,
                    15: 0.0422 }



if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument("-t", "--scoretype", required=True, choices=["overlap","nonoverlap"], help="")
    parser.add_argument('-p', "--selected_frags_path", required=True, default="", help='')
    parser.add_argument("-s", "--scorefile", help="")
    parser.add_argument("-w", "--closab_weight", action="store_true", default=False, help="when it says it's closable, reweight it")
    parser.add_argument("-nw", "--nonclosab_weight", type=float, default=10, help="when it says not closable, reweight it")
    args = parser.parse_args()

    # try to find indexing file from the folder, if not, create one
    if not args.scorefile:
        rsn = basename( getcwd() )
        args.scorefile = "pos_%s_score.sc" % rsn
        assert exists( args.scorefile )

    indexing = ScoreTable( args.selected_frags_path )

    ## read scorefile
    if args.closab_weight:
        outfn = args.scorefile+".clean.weighted.idx"
    else:
        outfn = args.scorefile+".clean.idx"

    lines = ""


    if args.scoretype == "overlap":
        lines = indexing.raw_overlap_score_cleaner( args.scorefile )

    elif args.scoretype == "nonoverlap":
        if args.closab_weight:
            lines = indexing.raw_nonoverlap_score_cleaner( args.scorefile, gap_weights_dict, args.nonclosab_weight )
        else:
            lines = indexing.raw_nonoverlap_score_cleaner( args.scorefile )

    if not lines:
        stderr.write("ERROR: nothing being done!\n"); exit()

    outfile = open( outfn, "w" )
    outfile.write( lines )
    outfile.close()

