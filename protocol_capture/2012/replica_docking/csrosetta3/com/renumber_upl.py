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
import string
#toolbox headers
import library

parser = ExampleArgumentParser(prog=basename(__file__),
                               description="Trim and renumber upl (NOE upper limit) file for given residue range",
                               examples=['%(prog)s final.upl trimmed.upl -s 5 -e 20'])
parser.add_argument("infile", help="chemical shift file");
parser.add_argument("outfile", help="chemical shift file",nargs="?",default="stdout");
parser.add_argument("-s","--start",dest="start",default="1",type=int, help="starting residue");
parser.add_argument("-e","--end",dest="end",default="0",type=int,help="ending residue");
parser.add_argument("-fasta",dest="fasta",help="figure out trimming from given sequence");
library.add_standard_args( parser )

args = parser.parse_args()

#input:
lines = open( args.infile,'r').readlines();

#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello( __file__ )

fasta=None
if args.fasta:
    fasta=library.read_fasta(args.fasta)

sequence="";
end=args.end
if end<=0:
	end=1000000
start=args.start
format="%5d %5s %5s   %5d %5s %5s  %8.3f     %s\n"
try:
	if fasta:
		upl_fasta=library.upl2fasta( lines )
		offset=library.find_fasta_offset( upl_fasta, fasta, verbose )
		start=offset+1
		end=start+len(fasta)-1

	for line in lines:
		#print line
		tags=string.split(line);
		if len(tags)<8: continue
		if int(tags[0])>=start and int(tags[0])<=end and int(tags[3])<=end and int(tags[3])>=start:
			outfile.write(format%(int(tags[0])-start+1,tags[1],tags[2],int(tags[3])-start+1,tags[4],tags[5],float(tags[6])," ".join(tags[7:])))

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
