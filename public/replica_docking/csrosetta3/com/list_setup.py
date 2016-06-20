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
import os
import traceback
import argparse
import sys

### toolbox library
try:
	import automatic_setup
	from library import LibException, MissingInput
	import library

except ImportError as exc:
	traceback.print_exc(exc)
	print "\ncall 'source %s/init'"%path.dirname(__file__)
	print "before using the toolbox or put this call into your .bashrc"
	exit()



name=sys.argv[0]
if 'list_setup' in path.basename( name ):
	alias='setup_list.py'
else:
	alias='list_setup.py'


#setup common-application cmdline options
application = automatic_setup.SetupApplication( 'This tool lists all setups in the target-library. One can filter-out setups using options -method/-methods, '+
																'-label/-labels and -target/-targets. Use flag -detail to get the full information of the listed setups. '+
																'This application does not change any setup'+
																'(requires choice of -method) Alias: %s'%alias,__file__, sys.argv, method_integration=False )

#setup specific cmdline options
parser=application.parser
parser.add_argument("-targets","-target",dest='targets',nargs='*', help="you can filter the display to a subset of targets");
parser.add_argument("-labels",nargs='*', default=None)
parser.add_argument("-methods",nargs='*')
parser.add_argument("-details",help="show full details of all listed setups (requires option -method)", default=False, action='store_true')
parser.add_argument("-show", help="show details for selected options in table form (requires option -method)", nargs='*', default=None )
#check usage  -- no specific method options allowed
application.check_absence_method_options()

args=parser.parse_args()
if ( args.details or args.show ) and not application.method:
	print "Option -detail or -show only works when you choose a method with -method"
	exit()

#main program
try:
	#look in targetlib for subdirectories
	targets=os.listdir(args.target_prefix)
	if args.show:
		msg="%-15s %-15s "%("Target","Label") +" %-15s"*len(args.show)%(tuple(args.show))
		print msg
		print "-"*len(msg)
	for t in targets:
		#filter targets if args.targets is set -- this could be done by regexpression (glob ? )
		if args.targets and not t in args.targets: continue
		#remove things that are just a file and not a subdirectory
		if not path.isdir( args.target_prefix+"/"+t ): continue
#		#attempt to instantiate a TargetDir object
		target=automatic_setup.TargetDir( target=t, prefix=args.target_prefix )
		if not target.exists(): continue
		#at this point target is a valid target

		#look for subdirectories that are methods
		methods=os.listdir(target.dir())
		for m in methods:
			#filter: a) it must be in method_choices,
			#        b) cmdline method selection, single method
			#        c) cmdline method selection, method list
			if not m in application.method_choices: continue
			if application.method and m!=application.method.name: continue
			if args.methods and not m in args.methods: continue
			if args.method and not m in args.method: continue
			#okay, this is a valid method -- get the variant labels
			labels=os.listdir(target.dir()+"/"+m)
			for l in labels:
				#obtain clear method-instance
				application.fresh_method()
				#apply cmdline filters. a) list, b) single label
				if args.labels and not l in args.labels: continue
				if '-label' in sys.argv and l!=args.label: continue

				#instantiate Setup for this set (target, method, label )
				setup=automatic_setup.Setup(target,m,l)

				#prints a line like: Setup( target, method_label )
				if args.show:
					if application.method:
						set_args=application.load_setup(setup, verbose=False)
					vals = [ getattr(set_args,x) for x in args.show ]
					print "%-15s %-15s "%(t,setup.label) + " %-15s"*len(vals)%(tuple(vals))
				else:
					print setup

				#print details, (need to load them)
				if application.method and args.details:
					setup.print_file_list()
					application.load_setup(setup)

#handle exceptions
except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
#        print sys.exc_type
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
 #       print inst
