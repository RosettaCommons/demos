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
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

#application specific headers

#default headers
import argparse
from ExampleArgumentParser import ExampleArgumentParser
from os.path import basename
import traceback, sys

#toolbox headers
import library
import cs
parser = ExampleArgumentParser(prog=basename(__file__),
															 description="convert chemical shifts from .tab (TALOS) to .prot (CYANA) format",
															 examples=['%(prog)s cs.tab > resonances.prot','%(prog)s cs.tab resonances.prot' ])
parser.add_argument("infile", metavar='cs.tab', help="chemical shift file");
parser.add_argument("outfile", metavar='resonances.prot', help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-header", help="write also header into prot file", action='store_true', default=False );
parser.add_argument("-noheader", help="suppress printing of the header", action='store_false', dest='header', default=False )
library.add_standard_args( parser )
args = parser.parse_args()

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello( __file__ )

try:
#	sequence=None
#	if args.fasta:
#		sequence=library.read_fasta( args.fasta )
#	if args.seq:
#		sequence=library.read_aa3_sequence( args.seq )

	talos = cs.TalosCSFile()
	talos.read_file( args.infile )

	prot =cs. ProtCSFile()
	prot.from_table( talos.table, talos.sequence )
	prot.write( outfile, header=args.header)

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
