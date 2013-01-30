#!/usr/bin/env python2.7
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
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os.path import basename
import argparse
from ExampleArgumentParser import ExampleArgumentParser
import sys
from amino_acids import longer_names
import bmrb
### toolbox library
import library
from os.path import exists
import cs
import traceback
#############################

parser = ExampleArgumentParser(prog=basename(__file__), description="Extract chemical shifts from BMRB file and convert to .prot (AutoNOE) format",
examples=['%(prog)s 2jrm.bmrb 2jrm.prot',
					'%(prog)s 2jrm.brmb > 2jrm.prot'])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("bmrb", help="A bmrb file or the chemical shift section of an bmrb file");
parser.add_argument("outfile", metavar="prot", help="chemical shift file",nargs="?",default="stdout")
parser.add_argument("-header", help="write a header into the file", action='store_true', default=False )
library.add_standard_args( parser )

args = parser.parse_args()

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello(__file__)

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

try:
	bmrb_file, fasta=bmrb.read_cs_bmrb(args.bmrb,verbose)
	try:
		cs_data=bmrb_file.process_frame('assigned_chemical_shifts', bmrb.cs_loop_cols )
	except KeyError:
		raise library.InconsistentInput("File %s does not contain a frame of the category 'assigned_chemical_shifts'. Maybe not a proper BMRB file?"%args.bmrb)
	if verbose:
		if args.header: extra_msg='with header to '
		else: extra_msg=''
		print "writing chemical shifts %sfile in prot-format to %s..."%(extra_msg,args.outfile)

	cs_table=cs.NIH_table()
	cs_table.from_dict( cs_data )

	prot=cs.ProtCSFile()
	prot.from_table( cs_table, sequence=fasta )
	prot.write( outfile, header=args.header )


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)



