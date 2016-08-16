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

    assert isdir( args.input_dir ), args.input_dir
    ref_fragfiles = " ".join( glob("%s/after_rotation_frags.%s.%s.*pdb" %( args.input_dir, args.mers, rsn ) ) )
    assert not len(ref_fragfiles)==0

    fragfiles_list = []
    # saner
    j_rsds = range( 1, args.n_total_positions+1 )[rsn:]

    for pdb in glob("%s/after_rotation_frags.%s.*pdb" %( args.input_dir, args.mers )):
        frag_rsn = int( basename(pdb).split(".")[2] )

        # to prevent silly bug
        assert frag_rsn < args.n_total_positions, "ERROR: please reset the n_total_positions!\n"

        if frag_rsn in j_rsds:
            fragfiles_list.append( pdb )


    assert not len(fragfiles_list)==0
    fragfiles = " ".join( fragfiles_list )

    # get command line asesmbled
    print getcwd()
    config = denovo_model_building_util.read_config_file( args.config_file )

    PATH = config.get("path", "demo_dir")
    exe = PATH + "/rosetta/" + args.exe
    db  = PATH + "/rosetta/" + args.database


    #flags += "-mapfile %s\n" % args.mapfile
    flags = "-ref_fragfiles %s\n" % ref_fragfiles
    flags += "-fragfiles %s\n" % fragfiles
    flags += "-out:score pos_%s_score.sc\n" % rsn
    flags += "-clash_dist_threshold %s\n" % args.clash_dist
    flags += "-database %s\n" % db
    if args.native:
        flags  += "-native %s\n" % args.native

    arguments = slim_flags + "\n" + flags

    if args.noloophash:
        arguments += "-bound\n"
        arguments += "-bound_stringent\n"
    else:
        arguments += "-lh:db_path /work/kenjung/for/possu/3to25mer\n"

    if args.soften_clash_score:
        arguments += "-soften_clash_score\n"


    flags_fn = open("flags", "w")
    flags_fn.write( arguments + "-mute all\n" )
    flags_fn.close()
    #print "writing down args"
    system( "%s @flags" % exe )



if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--config_file', default="../denovo_model_building_scripts.cfg", help='')
    parser.add_argument('-d', '--input_dir', default="../../Step1_Place_fragments_into_density/candidate_fragment_placements/", help='')
    parser.add_argument('-c', '--clash_dist', default=2.0, type=float, help='')
    parser.add_argument('--noloophash', action="store_true", default=True, help='')
    parser.add_argument('--soften_clash_score', action="store_true", help='')
    #parser.add_argument('--mapfile', required=True, type=str, help='')
    parser.add_argument('--native', type=str, help='')
    parser.add_argument('-n', '--n_total_positions', default=1000, type=int, help='this just serves as a range toward the end, just a number that larger than total_rsds')
    parser.add_argument('--exe', default="cal_nonoverlap_score.static.linuxgccrelease", help='')
    parser.add_argument('--database', default="rosetta_database", help='')
    parser.add_argument('-m', '--mers', default=9, type=int, help='')
    args = parser.parse_args()

    run_jobs( args )
