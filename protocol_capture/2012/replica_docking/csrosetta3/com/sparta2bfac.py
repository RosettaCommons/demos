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
from library import square
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
                                 description="calculate chemical shift scores for given pdb-structures using Sparta+. "+
																 "Sparta+ generated files are stored under -pred_dir, and pdbs colored (via bfactors) according to "+
                                 "local fit of chemical shift scores can be generated using -bfac/nobfac with the name of the generated pdbs "+
																 "controlled via -prefix",
																 examples=[('%(prog)s S1.pdb S2.pdb S3.pdb -cs 2jrm.tab -bfac -prefix structures/bfac',
																					 'print chemshift scores to screen and generate bfac structures in the subdirectory structures/ with '+
																					 'file-names structures/bfac_S1.pdb structures/bfac_S2.pdb structures/bfac_S3.pdb'),
																					 ('%(prog)s S1.pdb -cs 2jrm.tab -selection 5-10+80-100',
																						'obtain chemical shift score for residues 5-10 and 80-100'),
																					 ('%(prog)s S1.pdb -cs 2jrm.tab -sigma 0.1',
																						'some SPARTA+ calculated sigma seem to be too low and lead to a domination of chemical-shift score '+
																						'and bfactor for the affected residues. Use the option -sigma to set the lower bound (default 0.3)')],
                                 add_help=True)


parser.add_argument("pdbs", nargs='*', help="pdb-file(s)");
parser.add_argument("-hm", help="store rmsd improvement due to ring-currents", action='store_true', default=False );
parser.add_argument("-prefix", help="prefix for output pdbs", default="bfac" )
parser.add_argument("-selection", help="residue selection", default="" )
parser.add_argument("-bfac", help="produce a bfactor structure [default]", action='store_true', default=True );
parser.add_argument("-nobfac", help="do not produce a bfactor structure", action='store_false', dest='bfac' );
parser.add_argument("-list_low_sigma", help="produce a list of residues that have a low sigma in one of the structures", action='store_true', default=False );
parser.add_argument("-sigma", type=float, help="lower bound for sigma (default=0.3)", default=0.3 );

parser.add_argument("-pred_dir", help="directory prefix for directories with SPARTA output", default='pred')
parser.add_argument("-overwrite", help="overwrite all existing SPARTA output", default=False, action='store_true')
parser.add_argument("-cs", dest='refCS', help='reference chemical shifts to compare with SPARTA predictions' )
library.add_standard_args( parser )
args = parser.parse_args()

verbose=True

####### program start
try:
	if verbose:
		library.hello( __file__ )
	# compute average chem-shift score for a cluster  -- check against single structure -- should give same
	# pre-set offset from multiple offsets...
	# put combination of atoms into a b-factor for these values: CS_DIFF, HM_SHIFT, (CS_DIFF+HM_SHIFT*0.6)^2-CS_DIFF^2
	low_sig_resi=None
	#setup the Sparta Session with refCS, scrap directory and overwrite
	spCalc=cs.SpartaCalculator( args.refCS, args.pred_dir, overwrite=args.overwrite )

	#run sparta for each pdb
	for pdb in args.pdbs:
		sparta=spCalc.get_sparta( pdb )

		#extract data columns
		cs_diff=sparta.get_slice( "CS_DIFF" )
		sigma=sparta.get_slice("SIGMA")
		hm_shift=sparta.get_slice("HM_SHIFT")

		#calculate cs_diff+hm_shift*0.6 to get the delta without the ring-currents
		cs_diff_noHM=cs.table_op( cs_diff, hm_shift, lambda x,y: x+y*0.6 )

		if args.hm:
			#obtain HM difference
			result=cs.table_op( cs_diff_noHM, cs_diff, lambda x,y: square(y)-square(x) )
		else:
		  #obtain score, either from deltaCS or from deltaCS_noHM  score=sum( deltaCS^2/sigma^2/4 )
			result=cs.table_op( cs_diff, sigma, lambda x,y: square(x)/square(max(y,args.sigma))/4 )
			result_HM=cs.table_op( cs_diff_noHM, sigma, lambda x,y: square(x)/square(max(y,args.sigma))/4 )

			#filter to restrict to selected atoms
			selection=cs.AtomSelection( args.selection )
	    sparta = cs.sum_table( result, selection )
	    sparta_HM = cs.sum_table( result_HM, selection )
	    print "model: %s  score: %5.3f noHM_score: %5.3f"%(pdb, sparta, sparta_HM )

			#write output structure if requested
			if args.bfac:
				output_name = basename( pdb );
				if '/' in args.prefix: library.mkdirp(args.prefix)
				if not '/'==args.prefix[-1]:
					output_name='_'+output_name
				cs.write_to_bfactor( pdb, "%s%s"%(args.prefix, output_name ), result, sum_residue=True, response=lambda x: sqrt(x) )

   	#create list of atoms for which predicted sigma is below threshold 'args.sigma' in one of the pdb structures
    #if requested by args
		#for each atom [key=(resi, name)] find the lowest sigma in any of the predictions (sigma dependes on structure)
		if args.list_low_sigma:
			if not low_sig_resi:
				low_sig_resi={}
				for key in sigma:
					low_sig_resi[key]=sigma[key]
			else:
				for key in sigma:
					if sigma[key]<low_sig_resi[key]:
						low_sig_resi[key]=sigma[key]


	if low_sig_resi:
		for key in low_sig_resi:
			if ( low_sig_resi[key] < args.sigma ):
				print key, low_sig_resi[key]


except Exception as inst:
    print traceback.print_exc( inst )
    print traceback.print_exception(sys.exc_type, sys.exc_value, None)
#       print inst

