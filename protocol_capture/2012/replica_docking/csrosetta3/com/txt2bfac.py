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
from warnings import *
#from PDBParser import PDBParser
from PDB import PDBIO
from library import square

#default headers
import argparse
from ExampleArgumentParser import ExampleArgumentParser
from os.path import basename
import traceback, sys

#toolbox headers
import library

parser = ExampleArgumentParser(prog=basename(__file__),
															 description="convert a column of numbers into bfac, if two columns first is assumed to be resiude number",
															 add_help=True)

parser.add_argument("-data", help="data which will go into bfactors");
parser.add_argument("pdbin", help="input pdb-file");
parser.add_argument("pdbout", help="output pdb-file");
library.add_standard_args( parser )
args = parser.parse_args()

#output:
verbose=1
if args.pdbout=="stdout":
	outfile=sys.stdout
	verbose=0
else:
	outfile=open(args.pdbout,'w');
	library.hello( __file__ )

verbose=True

####### program start
try:
	data=open(args.data,'r').readlines()
	two_column=( len(data[0].split())>1 )
	ct=1
	bfac={}
	for l in data:
		if two_column:
			t=l.split()
			bfac[ int( t[0]) ]=float(t[1])
		else:
			bfac[ ct ]=float(l)
		ct = ct + 1
		#end - for l in data:

	with catch_warnings(record=True) as w:
		parser=PDBParser()
		s=parser.get_structure( "t000_", args.pdbin )
		for chain in s[0]:
			for residue in chain:
				resid=residue.get_id()[1]
				for atom in residue:
					if resid in bfac:
						atom.set_bfactor( bfac[ resid ] )
					else:
						atom.set_bfactor( 0 )

		io=PDBIO()
    io.set_structure( s )
    io.save( args.o )

except Exception as inst:
    print traceback.print_exc( inst )
    print traceback.print_exception(sys.exc_type, sys.exc_value, None)
#       print inst

