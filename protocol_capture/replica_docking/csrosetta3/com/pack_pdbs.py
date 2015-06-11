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
from os.path import basename
from os import path
import traceback
import argparse
import sys
import os
import shutil
from glob import glob
import subprocess
### toolbox library
try:
	from library import LibException, MissingInput
	import library
	from ExampleArgumentParser import ExampleArgumentParser
except ImportError as exc:
	traceback.print_exc(exc)
	print "\ncall 'source %s/init'"%path.dirname(__file__)
	print "before using the toolbox or put this call into your .bashrc"
	exit()


parser = ExampleArgumentParser(prog=basename(__file__), description="pack individual structures into a single multi-model pdb-file",
examples=['%(prog)s S1.pdb S2.pdb > all.pdb',
          '%(prog)s -silent some.out > all.pdb',
					'%(prog)s -silent some.out -o all.pdb'])

#add special options that are only relevant for setup_run
parser.add_argument("pdbs", nargs='*', help="input pdb structure")
parser.add_argument("-o", help="output pdb structure", default="stdout")
parser.add_argument("-silent", help="read models from silent file", default=None)
library.add_standard_args( parser )

#parse cmd-line and check usage
args = parser.parse_args()

#start the main program
try:
	if len(args.pdbs)==0 and not args.silent:
		raise MissingInput("need to specify input data: either pdbs or -silent")


#output:
	verbose=1
	if args.o=="stdout":
		outfile=sys.stdout
		verbose=0
	else:
		outfile=open(args.o,'w');

####### program start
	if verbose:
		library.hello( __file__ )
	try:
		if args.silent:
			import random
			tmpdir_b='tmp_pack_pdbs_%s'%path.splitext(basename(args.silent))[0]
			attempts=0
			tmpdir=tmpdir_b+"_%d"%random.randint(1,1000000)
			while attempts<99 and path.exists( tmpdir ): #basic help to avoid clashing -- not a safe protection when multithreading though...
				tmpdir=tmpdir_b+"_%d"%random.randint
				attempts+=1
			if path.exists( tmpdir ): exit("Cannot create tmp-directory because it exists already")
			library.mkdirp( tmpdir )
			os.chdir( tmpdir )
			prog=library.rosetta_executable( 'ensemble_analysis' )
			if verbose:
				print 'Read silent file: %s...'%args.silent
			subprocess.call( '%s -in:file:silent ../%s -wRMSD 2 -out:pdb -out:level 200'%(prog,args.silent), shell=True )
			pdbs=glob('*pdb')
		else:
			pdbs=args.pdbs

		model_ct=0
		if verbose:
			print 'Pack models...'
		for pdb in pdbs:
			model_ct+=1
			outfile.write('MODEL     %4d\n'%model_ct)
			outfile.write('REMARK ROSETTA-TAG %s\n'%path.splitext(basename(pdb))[0])
			for line in open(pdb,'r'):
				outfile.write('%s'%line)
			outfile.write('ENDMDL\n')
	finally:
		if 'tmpdir' in locals():
			if verbose:
				print 'Remove tmp-directory: %s...'%tmpdir
			os.chdir( '..')
			shutil.rmtree(tmpdir)

except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
	exit(1)
