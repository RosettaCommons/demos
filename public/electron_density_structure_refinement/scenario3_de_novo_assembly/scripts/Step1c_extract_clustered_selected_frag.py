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
from os.path import basename, exists, isabs, getsize
from sys import exit, stderr, stdout
from glob import glob
import shutil

import denovo_model_building_util
from denovo_model_building_util import slim_flags

def get_clustered_results( clustered_results ):
    rsn = basename( getcwd() )
    candidate_placements = []
    with open( clustered_results, "r" ) as f:
        for l in f:
            if l.startswith("#"): continue
            ls = l.split()
            assert ls[0] == rsn
            candidate_placements.append( ls[-1].split(".pdb")[0] )

    return " ".join( candidate_placements )



def run_jobs( args ):
    denovo_model_building_util.touch("extract.lock")

    rsn = int( basename( getcwd() ) )
    tags = get_clustered_results( args.clustered_results )
    #print tags
    if not tags:
        stderr.write("ERROR: no tags available for this target\n")
        exit()

    #silent_fn = "pos_%s_silent_packed_c3_80percent_cut.out" % rsn
    #silent_fn = "pos_%s_silent_packed_c3_100percent_cut.out" % rsn
    matched_files = glob("pos_%s_%s*" %( rsn, args.silent_fn_pattern ))
    assert len(matched_files)==1, matched_files
    silent_fn = matched_files[0]

    config = denovo_model_building_util.read_config_file( args.config_file )
    PATH = config.get("path", "demo_dir")
    exe = PATH + "/rosetta/extract_pdbs.static.linuxgccrelease "
    db  = PATH + "/rosetta/rosetta_trunk_database"

    arguments = slim_flags + "-mute all -in:file:silent %s -in:file:silent_struct_type binary -in:file:tags %s -database %s" %( silent_fn, tags, db )
    cmd = exe + arguments

    if args.dry_run:
        print( cmd )
    else:
        #print cmd
        system( cmd )
        assert len( glob("after_rotation_frags.*.*.*.*.????.pdb" ) ) != 0
        #shutil.move("after_rotation_frags*pdb", args.outfrags_dir )
        system("mv after_rotation_frags*pdb %s" % args.outfrags_dir )



if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--config_file', default="../denovo_model_building_scripts.cfg", help='')
    parser.add_argument('-f', '--clustered_results', default="clustering_selected.sc", help='')
    parser.add_argument('-p', '--silent_fn_pattern', default="silent_5percent_cut", help='')
    parser.add_argument('-o', '--outfrags_dir', default="../../candidate_fragment_placements", help='')
    parser.add_argument('--dry_run', action="store_true", help='')
    args = parser.parse_args()

    assert exists( args.clustered_results ), "could not find %s, please run collect_clustering_results.sh" % args.clustered_results
    assert getsize( args.clustered_results ) > 0, "%s is empty" % args.clustered_results
    if not exists( args.outfrags_dir ):
        system("mkdir %s" % args.outfrags_dir )

    assert not exists("extract.lock"),"job has already existed. To rerun, remove extrack.lock"
    run_jobs( args )
