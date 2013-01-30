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

import string
#from sys import argv,stderr,exit
#import sys
from os.path import basename
from os import path
import sys
import traceback
import cs
import library
import shutil
import sets

from ExampleArgumentParser import ExampleArgumentParser

__dir__=path.split( __file__ )[0]
__topdir__=path.split(__dir__ )[0]

parser = ExampleArgumentParser(prog=basename(__file__),
															 fromfile_prefix_chars='@',
															 description="adjust chemical shift offset")


parser.add_argument("infile", help="talos file that should be adjusted");


group_individual=parser.add_argument_group('offset as individual options')
group_individual.add_argument("-N", help="adjustment for N", type=float,default=0)
groupCO=group_individual.add_mutually_exclusive_group()
groupCO.add_argument("-CO", dest='C', help="offset of CO", type=float,default=0)
groupCO.add_argument("-C", help="offset of CO", type=float,default=0)
group_individual.add_argument("-CA", help="offset of CA", type=float,default=0)
group_individual.add_argument("-CB", help="offset of CB", type=float,default=0)
groupHN=group_individual.add_mutually_exclusive_group()
groupHN.add_argument("-HN", help="offset of amide H", type=float,default=0)
groupHN.add_argument("-H", dest="HN", help="offset of amide H", type=float,default=0)
group_individual.add_argument("-HA", help="offset of HA", type=float,default=0)

parser.add_argument("-offsets", help="give offsets like this 'CB: 0.00ppm CA: -0.10ppm' or this 'CB 0.00 CA -0.10", default=None)

parser.add_argument('-noreverse',action='store_false',dest='reverse',
										help="do not reverse sign, use numbers directly as the adjustement",	default=True)
parser.add_argument('-keep_sign',action='store_false',dest='reverse',
										help="do not reverse sign, use numbers directly as the adjustement",	default=True)

group_outfile=parser.add_mutually_exclusive_group()
group_outfile.add_argument('-i',action='store_true',help='change in place (backup will be written with ~ as ending)', default=False)
group_outfile.add_argument("outfile", nargs='?', help="talos file with adjusted CS",default="stdout");
#parser.add_argument("-out", help="store averaged offset-corrected chemical shifts here" );
library.add_standard_args( parser )
args = parser.parse_args()


verbose=1
if args.outfile=="stdout":
    outfile=sys.stdout
    verbose=0
else:
    outfile=open(args.outfile,'w');

####### program start
try:
	if verbose:
		library.hello( __file__ )

 	from cs import TalosCSFile
	tab=TalosCSFile()
	tab.read_file( args.infile )

	sign=1
	if args.reverse:
		sign=-1

	untreated_atoms=sets.Set()
	offset_dict={ 'N':sign*args.N, 'CO':sign*args.C, 'C':sign*args.C, 'CA':sign*args.CA, 'CB':sign*args.CB, \
								'HN':sign*args.HN, 'H':sign*args.HN, 'HA':sign*args.HA,'HA1':sign*args.HA,'HA2':sign*args.HA,'HA3':sign*args.HA}
	ambiguities={'HN':['H','HN'], 'H':['HN','H'], 'HA':['HA','HA1','HA2','HA3'], 'CO':['C','CO'], 'C':['C','CO'] }
	if args.offsets:
		 #check that no direct offset is present
		 for (key,val) in offset_dict.iteritems():
			 if val != 0:
				 exit('cannot use -offsets and any of the individual offset options: -%s'%key)
		 tags=args.offsets.replace('ppm','').replace(':','').split()
		 if len(tags)/2*2!=len(tags):
			 exit('found uneven number of tags in -offsets string: require atom offset pairs')
		 for i in range(0,len(tags),2):
			 if not tags[i] in offset_dict:
				 exit('cannot recognize atom-name %s in -offsets string'%tags[i])
			 try:
				 to_change=ambiguities[tags[i]]
			 except KeyError:
				 to_change=[tags[i]]
			 for key in to_change:
				 offset_dict[key]=sign*float(tags[i+1])

	if verbose:
		print 'will apply the fullowing adjustments to the chemical shift table'
		for (key,val) in offset_dict.iteritems():
			print '%4s %6.2f'%(key,val)

	atom_col=tab.get_col('ATOMNAME')
	cs_col=tab.get_col('SHIFT')
	for entry in tab.table.itervalues():
		 try:
			 entry[cs_col]+=offset_dict[entry[atom_col]]
		 except KeyError:
			 untreated_atoms.add(entry[atom_col])
	if args.i:
		 backup=args.infile+'~'
		 ct=1
		 while path.exists(backup):
			 backup=args.infile+'.%d~'%ct
			 ct+=1
		 shutil.move(args.infile,backup)
		 tab.write_file( args.infile )
	else:
		 tab.write( outfile )

	if len(untreated_atoms) and verbose:
		print 'ignored atoms (no correction applied):'
		print " ".join(['%s'%atom for atom in untreated_atoms])

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)


#	cs_obs_av=None;
#	inv_N=len(args.pdbs);
#	for pdb in args.pdbs:
#		sparta=Sparta( pdb )
#		cs_obs=sparta.get_slice( "CS_OBS" )
#		if not cs_obs_av:
#			cs_obs_av=table_op( cs_obs, cs_obs, lambda x,y: x*inv_N )
#		else:
#			cs_obs_av=table_op( cs_obs_av, cs_obs, lambda x,y: x+y*inv_N )
# bla=sparta.pred
#	bla.write_file("manno")
