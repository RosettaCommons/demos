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
import Util
from sys import stdout, exit, stderr
from glob import glob
import numpy as np
from os import popen
from os.path import basename
from numpy.core.umath_tests import inner1d
import operator


def cal_rms( matA, matB ):
    L = matA.shape[0]
    if (L<1): return 0
    dx = matA - matB
    dx2 = inner1d(dx,dx)

    return np.sqrt( sum(dx2)/L )


def parse_results( args ):
    """ read through all the mc_results files
        return a dict[ rsn ] = { fragid: [ count, rmsd ], fragid2: [ count, rmsd ] ... } """

    dict = {}
    fragfiles = []
    count = 0
    for file in args.mc_results_fn:
        with open( file, "r" ) as f:
            #print file
            for l in f:
                ls = l.split()
                try:
                    rsn = int( ls[2].replace(",","") )
                except ValueError: # if total_score < 0
                    break

                rmsd = float( ls[4].replace(",","") )
                fragid = ls[7]

                if fragid == "null": continue

                fragfiles.append( fragid )
                fragid = Util.extract_fragid( ls[7] )

                if dict.has_key( rsn ):
                    if dict[ rsn ].has_key( fragid ):
                        dict[ rsn ][ fragid ][0] += 1
                    else:
                        dict[ rsn ][ fragid ] = [ 1, rmsd ]
                else:
                    dict[ rsn ] = { fragid : [ 1, rmsd ] }
        count += 1

        fragfiles = list( set( fragfiles ) )

    #print fragfiles
    return dict, fragfiles, count # total number of lines


def get_CA_matx( pdb ):
    '''
    get a dict = { residue_num:{ atom:xyz }}
    '''
    xyzDict = {}
    nres = 0
    # parse pdb xyz coordinates
    #print "get_CA_matx", pdb
    with open( pdb, "r" ) as f:
        for l in f:
            atom = l[12:16].strip()
            if atom != "CA":
                continue
            res   = int( l[22:26])
            xyz   = map(float, [ l[30:38], l[38:46], l[46:54] ])
            xyzDict[ res ] = xyz
            nres += 1


    # convert them into array for fast rmsd calculation
    matx = np.zeros([ nres, 3 ])
    rsds = sorted( xyzDict.keys() )
    offset = rsds[0]-1 # to prevent some pdbs doens't start with 1, but with their original numbering
    for res in rsds:
        #print res, res-offset-1
        #print xyzDict[ res ]
        matx[ res-offset-1 ] = xyzDict[ res ]

    return matx


def get_frag_xyzMatx( fragfiles, frags_dir ):
    """ return a dict[ fragid ] = len*3 xyz numpy array """
    dict = {}

    for frag in fragfiles:
        #print
        frag = frag.split("????.pdb")[0]
        #print frag
        frag_fn = glob( frags_dir+"/%s*pdb" % frag )
        #print frag_fn
        assert ( len( frag_fn ) == 1 ), frag_fn
        fragid = Util.extract_fragid( basename( frag_fn[0] ) )
        dict[ fragid ] = get_CA_matx( frag_fn[0] )

    return dict


def get_highest_count_fragid( dict ):
    """ return the fragid that has the most count in a position """
    fragid, [ count, rmsd ] = sorted( dict.iteritems(), key=operator.itemgetter(1) )[-1]
    return fragid, count, rmsd


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument("-f", "--mc_results_fn", nargs="+", required=True )
    parser.add_argument("-d", "--frags_dir", required=True )
    parser.add_argument("--confidence_cut", default=1.00, type=float, help="serves as a threshold to decide what percentage of total count that I would like retain" )
    parser.add_argument("--rmsd_threshold_to_add", default=3.0, type=float, help="rmsd threshold to add count to the representative_fragid" )
    args = parser.parse_args()

    # parse the mc results files
    mc_selected_frags_dict, fragfiles, total_count = parse_results( args )
    # get fragXYZ coordinates for each selected frags
    frag_matx_dict = get_frag_xyzMatx( fragfiles, args.frags_dir )

    # the count that I accept as "common"
    accepted_count = total_count*args.confidence_cut
    assert total_count >= accepted_count
    stdout.write("#total_count: %s; accepted_count: %s\n" %( total_count, accepted_count ) )

    for rsn in sorted( mc_selected_frags_dict.keys() ):
        # start a new outline
        outline = "%3s" % str(rsn)

        # reprt_fragid, the fragid that has the highest count
        reprt_fragid, orig_count, reprt_rmsd = get_highest_count_fragid( mc_selected_frags_dict[ rsn ] )
        new_count = orig_count

        #print mc_selected_frags_dict[rsn],
        for fragid in mc_selected_frags_dict[ rsn ].keys():
            #print fragid,
            if fragid != reprt_fragid:
                rms = cal_rms( frag_matx_dict[ reprt_fragid ], frag_matx_dict[ fragid ] )
                #stdout.write(" %s %s %s" %( reprt_fragid, fragid, rms ) )
                if rms <= args.rmsd_threshold_to_add:
                    new_count += mc_selected_frags_dict[ rsn ][ fragid ][0]

        outline += "  count: %3s" % new_count
        outline += "  orig_count: %3s" % orig_count
        outline += "  rmsd: %5.4s" %( reprt_rmsd )
        outline += "  fragfile: after_rotation_frags.%s.%s.%s.%s.????.pdb" % reprt_fragid
        outline += "\n"
        if new_count >= accepted_count:
            stdout.write( outline )

