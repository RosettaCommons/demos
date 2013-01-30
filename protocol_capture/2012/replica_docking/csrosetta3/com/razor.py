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

#from sys import argv,stderr,exit
#import sys
from os.path import exists
from os.path import basename
from os import path
from glob import glob
import traceback
import argparse
import sys
import os

from ExampleArgumentParser import ExampleArgumentParser
### toolbox library
try:
	import automatic_setup
	from library import LibException, MissingInput
	import library
except ImportError as exc:
	traceback.print_exc(exc)
	print "\ncall 'source %s/init'"%path.dirname(__file__)
	print "before using the toolbox or put this call into your .bashrc"
	exit(1)


parser = ExampleArgumentParser(prog=basename(__file__),
                               description="Remove residues from PDB without renumbering",
                               examples=['%(prog)s input.pdb -remove 1+2-5 > output.pdb',
																				 '%(prog)s input.pdb output.pdb -keep 6-50' ])


parser.add_argument("infile", metavar='in.pdb', help="input pdb file");
parser.add_argument("outfile", metavar='out.pdb',help="output pdb file",nargs="?",default="stdout");
parser.add_argument("-remove",help="residue range to remove (1-5+75-)", default=None);
parser.add_argument("-keep", help="residue range to keep (6-75)", default=None);
parser.add_argument('-rigid', help='read range from ROSETTA .rigid file', default=None);
library.add_standard_args( parser )

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def find_max_res( lines ):
	max_res=0
	for line in lines:
    line_edit = line
    if line[0:3] == 'TER':
        continue

    if line_edit[0:4] == 'ATOM' or line_edit[0:6] == 'HETATM':
			if not (line[16]==' ' or line[16]=='A'): continue
			resnum = int(line_edit[23:26])
			if resnum>max_res:
				max_res=resnum
	return max_res

def parse_ranges( args ):
	removed=[]
	if args.remove:
		if args.keep or args.rigid:
			raise library.InconsistentInput('Cannot choose simultaneously more than one of the options -remove, -keep or -rigid')
		ranges=args.remove.split('+')
		for r in ranges:
			start_end=r.split('-')
			try:
				start=int(start_end[0])
			except:
				start=1
			try:
				end=int(start_end[1])
			except:
				end=max_res
			for i in range(start,end+1):
				removed.append(i)
	if args.rigid:
		if args.keep or args.remove:
			raise library.InconsistentInput('Cannot choose simultaneously more than one of the options -remove, -keep or -rigid')
		rigids=open(args.rigid,'r').readlines()
		removed=range(1,max_res+1)
		for r in rigids:
			tags=r.split()
			if len(tags)==0: continue
			if tags[0]!='RIGID':
				raise library.InconsistentInput('Incorrect file format in %s. Expected RIGID at beginning of line'%args.rigid)
			for i in range(int(tags[1]),int(tags[2])):
				try:
					removed.remove(i)
				except:
					pass
	if args.keep:
		removed=range(1,max_res+1)
		if args.rigid or args.remove:
			raise library.InconsistentInput('Cannot choose simultaneously more than one of the options -remove, -keep or -rigid')
		ranges=args.keep.split('+')
		for r in ranges:
			start_end=r.split('-')
			try:
				start=int(start_end[0])
			except:
				start=1
			try:
				end=int(start_end[1])
			except:
				end=max_res
			for i in range(start,end+1):
				try:
					removed.remove(i)
				except:
					pass
	return removed



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
	pdblines=open(args.infile,'r').readlines()
	#parse lines to figure out max-res
	max_res=find_max_res( pdblines )

	removed=parse_ranges( args )
	print removed

	for line in pdblines:
    if line[0:3] == 'TER':
        continue

    if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
#			if not (line[16]==' ' or line[16]=='A'): continue
			resnum = line[23:26]
			if not int(resnum) in removed:
				outfile.write(line)

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

