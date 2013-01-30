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
import argparse
import sys
import library
import cs
from ExampleArgumentParser import ExampleArgumentParser
from os.path import basename
import traceback
#############################
#if len(argv) <=1:
#    Help()

#file = argv[1]
parser =ExampleArgumentParser(prog=basename(__file__),
										description="find rigid/flexible residues according to S2 RCI index in pred.tab of TALOS+ output",
										add_help=True)
parser.add_argument("infile",help="pred.tab")
parser.add_argument("outfile", help="rigid file",nargs="?",default="stdout");
parser.add_argument("-threshold",help="residues below this S2 value are considered flexible", type=float, default=0.7)
parser.add_argument("-cutoff",dest='threshold',help="residues below this S2 value are considered flexible", type=float, default=0.7)
parser.add_argument('-smooth', help='avoid tagging <= N consecutive residues as flexible', default=2, type=int )
parser.add_argument('-flexible', action='store_true',help='write flexible residues instead of rigid', default=False )
library.add_standard_args( parser )
args = parser.parse_args()

#take list of ranges: make string with RIGID XXX XXX
def range2rigid(ranges,tag='RIGID'):
	s=''
	for r in ranges:
		s+='%s %5d %5d\n'%(tag,r[0],r[1])
	return s

def invert_ranges(ranges,max_res):
	start=1
	inverted=[]
	for r in ranges:
		if r[0]-1>start:
			inverted.append((start,r[0]-1))
		start=r[1]+1
	if start<max_res:
		inverted.append((start,max_res))
	return inverted
#output:
verbose=1
if args.outfile=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.outfile,'w');
	library.hello( __file__ )

try:
	pred=cs.NIH_table()
	pred.read_file(args.infile)
	tab_s2=pred.get_slice('S2')
#	tab_count=pred.get_slice('COUNT')

	min_res=sorted(tab_s2.keys())[0]
	max_res=sorted(tab_s2.keys())[-1]
	#flex=[ resi for ((resi,s2),count) in zip(tab_s2.iteritems(),tab_count.itervalues()) if s2<args.threshold or count==0 ]
	#I looked at a bunch of files, s2==1.0 is always in flexible regions. ~0.7
	flex=range(1,min_res)+[ resi for (resi,s2) in tab_s2.iteritems() if s2<args.threshold or s2==1.0 ]

	try:
		start=flex[0]
		last=flex[0]
		ranges=[]
		for resi in flex:
			if resi>last+args.smooth+1:
				ranges.append((start,last))
				start=resi
			last=resi
		ranges.append((start,last))
		ranges=[r for r in ranges if r[1]-r[0] >= args.smooth ]
		smoothed_flex=sum([range(r[0],r[1]+1) for r in ranges ],[])
		if verbose:
			print 'flexible residues:      '+" ".join([ '%d'%f for f in flex])
			print 'smoothed flex residues: '+" ".join([ '%d'%f for f in smoothed_flex])
	except IndexError:
		pass

	if args.flexible:
		#ROSETTA uses LOOP rather than FLEX for flexible residues
		outfile.write(range2rigid(ranges,'LOOP'))
	else:
		outfile.write(range2rigid(invert_ranges(ranges,max_res)))

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
