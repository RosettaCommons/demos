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
from math import sqrt

#default headers
import argparse
from ExampleArgumentParser import ExampleArgumentParser
from os.path import basename
import traceback, sys

#toolbox headers
import cs
import library


parser = ExampleArgumentParser(prog=basename(__file__),
															 fromfile_prefix_chars='@',
															 description="Compute difference between two sets of chemical shifts (TALOS format) and write as bfactor to pdb files. "+
															 'Differences are summed up on per residue basis, and are weighted according to Evenas J, et al. 2001 JMB with the following weights '+
															 '"N":0.154,"H":1.0,"C":0.341,"CA":0.276, "CB":0.276.',
															 examples=[('%(prog)s -cs cs_apo.tab cs_holo.tab S1.pdb S2.pdb',
																					'store CSP perturbation as bfactor (per residue) in provided pdb files'),
																				 ],
															 add_help=True)

parser.add_argument("pdbs", nargs='*', help="pdb-file(s)");
parser.add_argument("-cs", nargs=2, help="compute difference between chemical shifts and store as bfac",required=True );
parser.add_argument("-prefix", help="prefix for output pdbs", default="bfac" )
parser.add_argument("-selection", help="residue selection, unselected residues will be set to 0", default=None )
parser.add_argument("-txt", help="dump CSP to txt-file with columns 'resid csp'", default=None )
library.add_standard_args( parser )
args = parser.parse_args()

verbose=True

####### program start
atom_weights = {"N":0.154,"H":1.0,"C":0.341,"CA":0.276, "CB":0.276 } #from Evenas J, et al. 2001 JMB
try:
	if verbose:
		library.hello( __file__ )

	#read chemical shift files
	csA=cs.NIH_table().read_file(args.cs[0]).get_slice("SHIFT")
	csB=cs.NIH_table().read_file(args.cs[1]).get_slice("SHIFT")

	#obtain residue selection
	selection=None
	if args.selection:
		sel=cs.AtomSelection( args.selection )

	#calculate difference for selected residues
	result = cs.table_op( csA, csB, lambda x,y : x-y, selection=sel )

	#store bfactors in all provided pdbs
	for pdb in args.pdbs:

		#output filename
		outfile="%s_%s"%(args.prefix,basename( pdb ))
		print "write deltaCS as bfactors to %s..."%outfile

		#compute square of difference (func), sum with given weights over all atoms of residue (sum_reside) and take the sqrt (response)
		delta_per_res=cs.write_to_bfactor( pdb, outfile, result, sum_residue=True, atoms=atom_weights, func=lambda x: x*x, response=lambda x: sqrt(x) )

	#if dump-is requested, write values to file
	if args.txt and delta_per_res:

		#output filename
		outfile=open(args.txt,'w')
		print "write deltaCS as plain text to file %s"%args.txt

		#dump
		for r,delta in delta_per_res.iteritems():
			outfile.write("%8d %8.3f\n"%(r,delta))

	#dump requested, but values not calculated
	elif args.txt:
		print "Structure required to compute deltaCS..."


except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

