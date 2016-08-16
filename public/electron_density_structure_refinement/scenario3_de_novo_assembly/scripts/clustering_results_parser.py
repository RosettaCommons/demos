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

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-r", "--results_file", required=True, help="")
    parser.add_argument("-s", "--silent_file", required=True, help="this is for getting header to figure out column")
    parser.add_argument("-c", "--cluster_size_cut", type=int, default=2, help="")
    parser.add_argument("-b", "--bucket_size", type=int, default=25, help="")
    parser.add_argument("--random_fill", action="store_true", default=False, help="")
    args = parser.parse_args()

    score_Dict = {}
    fields_Dict = denovo_model_building_util.read_silent_header( args.silent_file )
    for line in popen("grep ^SCORE: %s" % args.silent_file, "r" ).readlines():
        if "description" in line: continue
        line_ed = line.split()
        score        = float( line_ed[ fields_Dict["elec_dens_whole_structure_allatom"] ] )
        super_rmsd   = float( line_ed[ fields_Dict["after_super_rmsd"] ] )
        nosuper_rmsd = float( line_ed[ fields_Dict["after_no_super_rmsd"] ] )
        tag          = line_ed[ fields_Dict["description"] ]

        score_Dict[ tag ] = ( score, nosuper_rmsd, super_rmsd )

    #pprint( score_Dict )

    Dict = {}
    for line in open( args.results_file, "r" ).readlines():
        if line.startswith("test_cluster:"):
            line_ed    = line.split()
            mothership = line_ed[7]
            jet        = line_ed[4]
            if mothership not in Dict.keys():
                Dict[ mothership ] = { jet:score_Dict[ jet ] }
            else:
                Dict[ mothership ][ jet ] = score_Dict[ jet ]
    #pprint( Dict )
    for key in Dict.keys():
        print "#\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%5s\t%s:\t%s" %( min([ Dict[key][jet][1] for jet in Dict[key].keys() ] ),
                                                               max([ Dict[key][jet][1] for jet in Dict[key].keys() ] ),
                                                               min([ Dict[key][jet][0] for jet in Dict[key].keys() ] ),
                                                               max([ Dict[key][jet][0] for jet in Dict[key].keys() ] ),
                                                               len(Dict[key]),
                                                               key,
                                                               Dict[key] )

    bucket = []
    index = 1
    # filling the bucket by energy order or random?
    if args.random_fill:
        cluster_centers = Dict.keys()
    else:
        print "# cluster centers sorted by energy"
        cluster_centers = [ tuple[0] for tuple in sorted( [ ( ct, min([ Dict[ct][jet][0] for jet in Dict[ct].keys() ] ) ) for ct in Dict.keys() ], key=lambda value: value[1]) ]

    print "# cluster centers: ", cluster_centers
    assert len(cluster_centers) > 0 # silly empty file

    while len( bucket ) <= args.bucket_size:
        for ms in cluster_centers: # is the keys random enough? or should I sorted by energy?
            if len( Dict[ ms ] ) <= args.cluster_size_cut: continue
            #print ms
            sorted_Tuple = sorted( Dict[ ms ].items(), key=lambda value: value[1][0], reverse=1 )
            try:
                picked = sorted_Tuple[-index]
                bucket.append( picked )
                print "#", ms, picked
            except:
                pass
            if len( bucket ) >= args.bucket_size: break
        index+= 1
        if len( bucket ) >= args.bucket_size: break

    #pprint( bucket )
    for wow in bucket:
        print "%3s\t%6.3f\t%6.3f\t%6.3f\t%s" %( wow[0].split(".")[2], wow[1][0], wow[1][1], wow[1][2], wow[0]+".pdb" )



