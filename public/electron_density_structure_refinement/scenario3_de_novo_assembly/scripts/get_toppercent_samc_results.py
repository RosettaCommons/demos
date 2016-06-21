#!/usr/bin/env python2.7
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) All the files in this directory and sub-directories are part of the Rosetta software
# (c) suite and are made available under license.  The Rosetta software is developed by the
# (c) contributing members of the Rosetta Commons. For more information, see
# (c) http://www.rosettacommons.org. Questions about this can be addressed to University of
# (c) Washington UW TechTransfer, email: license@u.washington.edu.
#
#  @author Ray Yu-Ruei Wang, wangyr@u.washington.edu
#
import argparse
from sys import exit, stderr, stdout, path
from os import popen, system, path

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Input args
    parser.add_argument('-p', '--percentage', default=5, type=int, help='')
    parser.add_argument('-r', '--results_fn', default="SAMC_results.txt", help='')
    parser.add_argument('-c', '--count', type=int, help='')
    args = parser.parse_args()

    # get SCORE: from all mc_sampling.out
    cmd = "grep SCORE: ./*/mc_sampling.out > %s" % args.results_fn
    system(cmd)
    assert path.getsize( args.results_fn ) > 0, "make sure you have something in your mc_sampling.out\n"

    if not args.count:
        n_total_models = int( popen("wc -l %s" % args.results_fn ).readline().strip().split()[0] )
        count = n_total_models*args.percentage/100
    else:
        count = args.count

    stderr.write("%s score lines (models) are extracted\n" %count )
    outfile = str(args.percentage)+"percent_selected_models.txt"
    system("sort -nk 3 %s > tmp" % args.results_fn)
    system("head -n %s tmp > %s" %( count, outfile ))
    system("rm -f tmp")
