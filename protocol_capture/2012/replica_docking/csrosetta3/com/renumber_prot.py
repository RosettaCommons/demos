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
import traceback

### toolbox library
import library
import fasta

parser = argparse.ArgumentParser(prog=basename(__file__), description="renumber residues in chemical shift file")
parser.add_argument("infile", help="chemical shift file");
parser.add_argument("outfile", help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-s","--start",dest="start",default=None,type=int, help="starting residue");
parser.add_argument("-e","--end",dest="end",default=None,type=int,help="ending residue");
parser.add_argument('-correct_fasta', default=None, help='if there is no sequence information in the .prot file, you can supply sequence with this flag')
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
	target_fasta=0

	sequence="";
	start=args.start
	end=args.end

	from cs import ProtCSFile
	tab=ProtCSFile()
	tab.read_file( args.infile )
	sequence=tab.sequence
	if not sequence and args.correct_fasta:
		sequence=fasta.read_fasta(args.correct_fasta)
		tab.set_sequence(sequence)
	elif sequence and args.correct_fasta:
		sys.stderr('WARNING: overwriting sequence in .prot file with input from -correct_fasta is this really intended?\n')
	if args.fasta:
		if start or end:
			exit('cannot choose -fasta together with -start and -end for trimming')
		if sequence:
			target_fasta=fasta.read_fasta(args.fasta)
			start=-fasta.find_fasta_offset(target_fasta,sequence,verbose)+1
			end=start+len(target_fasta)-1;
		else:
			exit('WARNING: cannot use fasta to trim since there is no sequence information in the .prot file')

	if args.rigid:
		if start or end: exit('cannot choose -fasta together with -start and -end for trimming')
		start,end = library.read_rigid_file( args.rigid )

	if not start:
		start=1
	if not end:
		end=0

	if verbose>0 and end: print 'Will trim from %d to %d'%(start,end)
	if verbose>0 and not end: print 'Will trom from %d to end'%(start,end)
	if sequence:
		sequence, end=fasta.cut_sequence(sequence,start,end,verbose)

	tab.renumber( start, end )

	tab.write( outfile )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)


