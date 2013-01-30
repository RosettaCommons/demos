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

parser = ExampleArgumentParser(prog=basename(__file__),
															 description="Convert upper limit restraint file (upl) into rosetta format. "+
															 "Cyana output final.upl usually contains ranking of restraints with SUP. "+
															 "it is advisable to separate the certain assignments SUP=1 from the less certain ones. "+
															 "Use option -low_qf for this purpose. If the target sequence is supplied "+
															 "(either in .fasta format or in .seq format) it will be used to automatically work out the offset "+
															 "from upl-numbering to the target sequence.",
															 examples=[('%(prog)s final.upl highQF.cst -low_qf lowQF.cst -fasta t000_.fasta',
																					'parse final.upl and write restraints with SUP=1 into highQF.cst and restraints with SUP<1'+
																					'into lowQF.cst. Use the target sequence in t000_.fasta to figure out trimming and offset'),
																				 ('%(prog)s final.upl all.cst','convert all restraints from final.upl to rosetta format')]
															 )
parser.add_argument("infile", metavar='upl', help="CYANA upl file");
parser.add_argument("outfile", metavar='cst', help="ROSETTA cst file",nargs="?",default="stdout");
parser.add_argument("-low_qf", help="write all restraints with SUP<1 into this file instead of main-output file", default=None )
parser.add_argument("-seq", help="sequence file to find offset and trimming");
parser.add_argument("-fasta", help="fasta file to find offset and trimming");
parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )
parser.add_argument("-min_sep", help="set N such that only restraints with i to i+N or more are written to output (default 4)", type=int, default=4)
parser.add_argument("-pad", help="pad upper distance bound by x (default 0.15)", type=float, default=0.15)

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
	sequence=None
	if args.fasta:
		sequence=library.read_fasta( args.fasta )
	if args.seq:
		sequence=library.read_aa3_sequence( args.seq)

	library.upl2mini( args.infile, outfile, fasta=sequence, QFall_file=args.low_qf, sep=args.min_sep, pad=args.pad, verbose=verbose )

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

