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
from os import popen, system
from os.path import basename, exists
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
import random

def bin( args ):
    dic = {}
    with open( args.results_file, "r" ) as f:
        #print file
        for l in f:
            ls = l.split()
            model = ls[0].split(":")[1]
            score = float( ls[2].replace(",","") )
            correct_rate = float( ls[4].replace(",","") )
            if score not in dic.keys():
                dic[ score ] = { model : correct_rate }
            else:
                dic[ score ][ model ] = correct_rate

    return dic



if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument("-f", "--results_file", required=True, help="")
    parser.add_argument("-o", "--output_folder", action="store_true", help="")
    args = parser.parse_args()
    Dict = bin( args )

    if args.output_folder:
        output_folder = args.results_file.split(".txt")[0] +"/"
        if not exists( output_folder ):
            system("mkdir %s" % output_folder )
    else:
        output_folder = "./"

    for score in Dict.keys():
        #print x, Dict[x][ random.sample( Dict[x] ,1) ]
        selected_model = random.sample( Dict[score].keys(), 1)[0]
        run = selected_model.split("_")[1]
        #print x, run, selected_model, ">"
        print "grep '%s:' %s/mc_sampling.out > %s%s.txt" %( selected_model, run, output_folder, selected_model )
        #print x, Dict[x]


