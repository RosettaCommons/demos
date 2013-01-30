#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
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
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen
from os import dup2
from os.path import exists
from operator import add
from math import sqrt
from os.path import basename
import argparse
import sys
from ExampleArgumentParser import ExampleArgumentParser
### toolbox library
import library
import traceback
import fasta


parser = ExampleArgumentParser(prog=basename(__file__),
															 description="renumber residues in chemical shift file",
                               examples=[('(%(prog)s in.tab out.tab -s 5 -e 76','keep residue 5 to 76, start counting in output at 1'),
																				 ('(%(prog)s in.tab -fasta target.trim.fasta > out.tab',
																					'figure out correct trimming from target.trim.fasta, trim output, start counting at 1 for output file')])

parser.add_argument("infile", help="chemical shift file");
parser.add_argument("outfile", help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-s","--start",dest="start",default="1",type=int, help="starting residue");
parser.add_argument("-e","--end",dest="end",default="0",type=int,help="ending residue");
mutex=parser.add_mutually_exclusive_group()
mutex.add_argument("-fasta",help="figure out trimming from given sequence");
mutex.add_argument("-rigid",help="use first and last rigid residue from .rigid file as written by pred2rigid",default=None)
library.add_standard_args( parser )
args = parser.parse_args()

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


#output:
verbose=1
if args.outfile=="stdout":
    outfile=sys.stdout
    verbose=0
else:
    outfile=open(args.outfile,'w');

####### program start
if verbose:
    library.hello( __file__ )

try:
	target=0

	sequence="";
	end=args.end

	from cs import TalosCSFile
	tab=TalosCSFile()
	tab.read_file( args.infile )
	sequence=tab.sequence
	start=args.start
	end=args.end

	if args.fasta:
    target=fasta.read_fasta(args.fasta)
		start=-fasta.find_fasta_offset(target,sequence)+1
		end=start+len(target)-1;

	if args.rigid:
		start,end = library.read_rigid_file( args.rigid )

	if verbose>0: print 'Will trim from %d to %d'%(start,end)
	if sequence:
    sequence, end=fasta.cut_sequence(sequence,start,end,verbose)

	tab.renumber( start, end )
	tab.write( outfile )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)


