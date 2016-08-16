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


def printer( rsd_list ):
    #print "# %s rsds: %s" %( len(rsd_list), " ".join( map( str, rsd_list ) ) )
    for rsd in rsd_list:
        print rsd

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("--pdb", required=True, help="")
    parser.add_argument("-f", "--fasta", required=True, help="")
    parser.add_argument("--assigned_rsds", action="store_true", default=False, help="")
    args = parser.parse_args()

    seq = ( denovo_model_building_util.fasta_file_reader( args.fasta ) )
    seq_rsds = range(1, len(seq)+1)

    xyz_dict, junk, pdbline_dict, junk  = denovo_model_building_util.create_xyzDict( args.pdb )
    assigned_rsds = xyz_dict.keys()
    unassigned_rsds = list( set(seq_rsds) - set(assigned_rsds) )

    if args.assigned_rsds:
        #print "# assigned_rsds"
        rsd_list = assigned_rsds
    else:
        #print "# unassigned_rsds"
        rsd_list = unassigned_rsds

    printer( rsd_list )
