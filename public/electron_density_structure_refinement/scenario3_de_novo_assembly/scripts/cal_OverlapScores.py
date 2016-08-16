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
from os import popen, system, getcwd
from os.path import basename, exists, isdir
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
from glob import glob

import denovo_model_building_util
from denovo_model_building_util import slim_flags



def run_jobs( args ):
    assert not exists("run.lock"), "job has already existed\n"
    denovo_model_building_util.touch("run.lock")

    rsn = int( basename( getcwd() ) )
    print rsn

    assert isdir( args.input_dir )
    ref_fragfiles = " ".join( glob("%s/after_rotation_frags.%s.%s.*pdb" %( args.input_dir, args.mers, rsn ) ) )
    assert not len(ref_fragfiles)==0

    fragfiles_list = []
    for pdb in glob("%s/after_rotation_frags.%s.*pdb" %( args.input_dir, args.mers )):
        this_rsn = int( basename(pdb).split(".")[2] )
        #if this_rsn in range( rsn, rsn+args.mers ) and this_rsn != rsn:
        if this_rsn in range( rsn, rsn+args.mers ) and this_rsn != rsn:
            fragfiles_list.append( pdb )

    assert not len(fragfiles_list)==0
    fragfiles = " ".join( fragfiles_list )

    print getcwd()
    config = denovo_model_building_util.read_config_file( args.config_file )

    PATH = config.get("path", "demo_dir")
    exe = PATH + "/rosetta/" + args.exe
    db  = PATH + "/rosetta/" + args.database

    #flags = " -native %s -mapfile %s -ref_fragfiles %s -fragfiles %s -out:score pos_%s_score.sc " %( args.native, args.mapfile, ref_fragfiles, fragfiles, rsn )
    flags = " -ref_fragfiles %s -fragfiles %s -out:score pos_%s_score.sc " %( ref_fragfiles, fragfiles, rsn )
    flags += " -overlap_dist_cutoff %s " % args.overlap_dist_cutoff
    flags += " -steepness_wt %s " % args.steepness_wt
    flags += " -seq_sep %s " % args.seq_sep
    flags += " -clash_dist %s " % args.clash_dist
    flags += " -clash_return_score %s " % args.clash_return_score
    flags += " -database %s " % db
    if args.native:
        flags += " -native %s " % args.native

    arguments = slim_flags + flags

    flags_fn = open("cmd.sh", "w")
    flags_fn.write( exe + arguments + " -mute all " )
    flags_fn.close()
    system( "sh cmd.sh" )



if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--config_file', default="../denovo_model_building_scripts.cfg", help='')
    parser.add_argument('-d', '--input_dir', default="../../Step1_Place_fragments_into_density/candidate_fragment_placements/", help='')
    parser.add_argument('--native',  help='')
    parser.add_argument('--mers', default=9, type=int, help='')

    parser.add_argument('--overlap_dist_cutoff', default=3.0, type=float, help='')
    parser.add_argument('--steepness_wt',        default=8.0, type=float, help='')
    parser.add_argument('--seq_sep',             default=5, type=int, help='')
    parser.add_argument('--clash_dist',          default=2.0, type=float, help='')
    parser.add_argument('--clash_return_score',  default=8.0, type=float, help='')

    parser.add_argument('-g', '--go_ahead', action="store_true", help='')
    parser.add_argument('--exe', default="cal_overlap_score.static.linuxgccrelease", help='')
    parser.add_argument('--database', default="rosetta_database", help='')
    args = parser.parse_args()

    run_jobs( args )
