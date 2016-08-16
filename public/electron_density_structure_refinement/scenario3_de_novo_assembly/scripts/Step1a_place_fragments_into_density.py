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
from os import system, getcwd
from os.path import basename, exists
from sys import exit, stderr

import denovo_model_building_util


def run_jobs( args ):
    ''' '''
    rsd_pos = basename( getcwd() )

    assert exists( args.flags_fn ), "flags does not existed\n"
    assert not exists("running.lock"), "job at %s has already launched, to relaunch please remove 'running.lock'\n" % rsd_pos

    denovo_model_building_util.touch("running.lock")

    # get all the necessary paths
    config = denovo_model_building_util.read_config_file( args.config_file )
    PATH = config.get("path", "demo_dir")
    exe = PATH + args.rosetta_dir + "/" + args.exe
    db  = PATH + args.rosetta_dir + "/" + args.database

    # get sequence position from the dir name to assemble arguments to run placements
    flags = " " + " ".join( [ line.strip() for line in open( args.flags_fn, "r" ).readlines() ] )
    flags += " -pos %s -out:file:scorefile pos_%s_score.sc -out:file:silent pos_%s_silent.out -database %s" %( rsd_pos, rsd_pos, rsd_pos, db )
    print( exe + flags + " -mute all" )
    system( exe + flags + " -mute all" )



if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--config_file', default="../denovo_model_building_scripts.cfg", help='')
    parser.add_argument('--flags_fn', default="../input_files/placement_flags", help='')
    parser.add_argument('--rosetta_dir', default="rosetta", help='')
    parser.add_argument('--exe', default="place_fragment_into_density.static.linuxgccrelease", help='')
    parser.add_argument('--database', default="rosetta_database", help='')
    opts = parser.parse_args()

    run_jobs( opts )
