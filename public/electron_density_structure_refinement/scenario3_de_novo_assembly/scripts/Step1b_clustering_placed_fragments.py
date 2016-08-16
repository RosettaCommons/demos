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


def extract_cmd( silentfile ):
    outfile = silentfile.split(".out")[0] + "_%spercent_cut.out" % str(opts.extract_percent).replace(".","")
    exe = PATH + "scripts/extract_lowscore_pdbs.py "
    args = "%s -p %s -so -ow -o %s --muteall" %( silentfile, opts.extract_percent, outfile )
    cmd = exe + args
    if not opts.dry_run:
        system( cmd )
    else:
        print cmd

    return outfile

def clustering_cmd( silentfile, outfile="clustering.log" ):
    exe = PATH + "rosetta/kcluster.static.linuxgccrelease "
    args = "-K_not_fit_xyz -K_level 1 -K_radius %s -in:file:silent %s -in:file:fullatom -in:file:silent_struct_type binary > %s" %( opts.cluster_radius, silentfile, outfile )
    cmd = exe + args

    if not opts.dry_run:
        system( cmd )
    else:
        print cmd

    return outfile

def parse_results( extract_out, clustering_out ):
    clustering_select = "clustering_selected.sc"
    exe = PATH + "scripts/clustering_results_parser.py "
    args = "-r %s -s %s -b %s > %s" %( clustering_out, extract_out, opts.bucket_size, clustering_select )
    cmd = exe + args

    if not opts.dry_run:
        system( cmd )
    else:
        print cmd


def run_jobs( args ):
    assert not exists( args.lockfile ), "%s has already existed\n" % args.lockfile

    if not args.dry_run:
        denovo_model_building_util.touch( args.lockfile )

    rsn = int( basename( getcwd() ) )

    if args.silentfn_pattern:
        silent_fn  = "pos_%s_silent_%s.out" %( rsn, args.silentfn_pattern )
    else:
        silent_fn  = "pos_%s_silent.out" % rsn

    print silent_fn
    assert exists( silent_fn )
    extract_out    = extract_cmd( silent_fn )
    clustering_out = clustering_cmd( extract_out )
    parse_results( extract_out, clustering_out )


if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument('--lockfile', default="clustering.lock", help='')
    parser.add_argument('--config_file', default="../denovo_model_building_scripts.cfg", help='')
    parser.add_argument("-r", "--cluster_radius", default=2.0, type=float, help='')
    parser.add_argument("-e", "--extract_percent", default=5, type=float, help='')
    parser.add_argument("-b", "--bucket_size", default=50, type=int, help='')
    parser.add_argument('-p', '--silentfn_pattern', default="", type=str, help='')
    parser.add_argument('-d', '--dry_run', action="store_true", help='')
    opts = parser.parse_args()

    config = denovo_model_building_util.read_config_file( opts.config_file )
    PATH = config.get("path", "demo_dir")

    run_jobs( opts )
