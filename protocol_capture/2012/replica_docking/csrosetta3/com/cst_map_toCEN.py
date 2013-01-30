#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
###
###
### This file is part of the CS-Rosetta Toolbox and is made available under
### GNU General Public License
### Copyright (C) 2011-2012 Oliver Lange
### email: oliver.lange@tum.de
### web: www.csrosetta.org
###
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <http://www.gnu.org/licenses/>.
###
###

import string
from glob import glob
from sys import argv,stderr,exit,stdout
from os import popen,system
from os.path import exists
from operator import add
from math import sqrt
import argparse
import traceback
#from pymol.vfont import plain
import sys
import library
import noe_tools
from noe_tools import cst_map_to_CEN_by_aa_type
import os

#############################
def Help():
    print '\n'
    print '-'*75
    print 'USAGE: %s fa_cst'%argv[0]
    print '\n will map all QX to CB and add 1.5 padding per methyl .. result printed to screen'
    print '\n\n\n\n'
    exit()

#if len(argv) <=1:
#    Help()


#file = argv[1]
parser =argparse.ArgumentParser(description="mapt the cst between H-H to the cst between C-C",
                                add_help=True)
parser.add_argument("-fasta",help="fasta file",required=True)
parser.add_argument("-library", help="library of CB-H dist",default=None )
#parser.add_argument("-c_type",help="which C atom is mapped",default="CB")
parser.add_argument("infile",help="(fullatom) constraint file")
parser.add_argument("outfile", help="(centroid) constraint file",nargs="?",default="stdout");
parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )

args = parser.parse_args()



#main program
try:
#output:
    verbose=1
    if args.outfile=="stdout":
        outfile=sys.stdout
        verbose=0
    else:
        outfile=open(args.outfile,'w');

    seq=library.read_fasta(args.fasta)

    #print seq
    cst_map_to_CEN_by_aa_type( args.infile, outfile, seq, mapping_library_file=args.library, verbose=verbose )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

