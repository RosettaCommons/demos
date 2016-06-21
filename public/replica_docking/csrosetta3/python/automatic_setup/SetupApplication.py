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

from glob import glob
#from sys import argv,stderr,exit
#import sys
from os.path import exists
from os.path import basename
import os
from os import path
import argparse
from ExampleArgumentParser import ExampleArgumentParser
import library
from Setup import Setup
from TargetDir import TargetDir

def extract_arg( argv, arg, default=None ):
	val=default
	if argv.count(arg):
		pos=argv.index(arg)
		val=argv[pos+1]
	return val


#common code for all setup applications:
# setup_target, setup_run, display_setup aka setup_display
class SetupApplication:
	def __init__( self, description, file, argv, method_integration=True, examples=None, aliases=None, setup_run=False ):
		self._target=None
		self.__file__ = file
		self.__dir__=path.split( self.__file__ )[0]
		self.__topdir__=path.split(self.__dir__ )[0]
		self.full_cmdline = argv

#output:
		verbose=1
		if verbose:
			library.hello(self.__file__ )


		self.flag_lib=extract_arg(argv,'-flag_lib',default=self.__topdir__+"/flag_library")

#now include the flag_lib/methods/XXX/options.py files
		self.method_choices=[]
		for methodf in glob( self.flag_lib+"/methods/*" ):
			if basename(methodf)[0]=="_": continue
			if 'TEMPLATE'==basename(methodf): continue
#			print "integrate method: %s from library"%basename(methodf)
			self.method_choices.append(basename(methodf))

			#setup main cmd-line parser
			self.parser = ExampleArgumentParser(prog=basename(self.__file__),
																					fromfile_prefix_chars='@',
																					description=description,
																					examples=examples,
																					aliases=aliases,
																					add_help=True)
		#load method
		self.method=None
		self.method_arg_group=None
		choice=extract_arg(argv, '-method' )
		has_help='-h' in self.full_cmdline or '--help' in self.full_cmdline
		if not has_help or method_integration:
			if choice in self.method_choices:
				method_code=self.flag_lib+"/methods/"+choice+"/options.py"
				method_name=choice
				method_path=self.flag_lib+"/methods/"+choice

				GLOBAL={"parser":self.parser,"method_name":method_name,"method_path":methodf,"flag_lib":self.flag_lib,"method":self.method}
				if setup_run:
					self.hidden_parser=argparse.ArgumentParser()
					group=self.hidden_parser.add_argument_group(title='method-options')
					run_group=self.parser.add_argument_group(title='run-options', description='extra options selected when generating a "RUN" with setup_run')
					GLOBAL['run_group']=run_group
					GLOBAL["group"]=group
				else:
					self.hidden_parser=None
					group=self.parser.add_argument_group(title='method-options', description='extra options for method')
					GLOBAL["group"]=group
				if exists( method_code ):
					exec open(method_code,'r') in GLOBAL
				else:
					print "CANNOT FIND METHOD CODE %s"%method_code
					exit()
				self.method = GLOBAL['method']
		#		print method_name, method_code
				self.method.name=method_name
				self.method_arg_group = None
				if 'group' in GLOBAL:
					self.method_arg_group = GLOBAL['group']

		#add other cmdline options
		self.add_common_options( self.__topdir__ )

	def add_common_options( self, top_dir  ):

		#check environment to set some defaults
		env_targetlib='CS3_BENCH_TARGETLIB'
		if env_targetlib in os.environ:
			target_prefix=os.environ[env_targetlib]
		else:
			target_prefix=None

		env_flaglib='CS3_BENCH_FLAGLIB'
		if env_flaglib in os.environ:
			flaglib=os.environ[env_flaglib]
		else:
			flaglib=top_dir+"/flag_library"

		#setup the cmdline options
		self.parser.add_argument("-target_dir", help="input-files are used from/put into this path, instead of directly into the run-directory");
		self.parser.add_argument("-target_prefix", help="target dirs are found here", default=target_prefix );
		self.parser.add_argument("-quiet", help="suppress output", action='store_true', default='false');
		self.parser.add_argument("-flag_lib", help="directory with meta files for generation of Rosetta cmd-line flag-files", default=flaglib );
		self.parser.add_argument("-label", help="give a special label to this setup -- reflected in file-name for method-options", default='standard' )
		self.parser.add_argument("-method", help="choose algorithm; use together with -h to see additional method-specific options", choices=self.method_choices );
		self.parser.add_argument("-remove", help="remove one or more options from the SETUP", nargs='*', action='append',default=[] )
		self.parser.add_argument("-patch", help="provide some python code that patches file-library content" );
		return self.parser


	#exit if no message has been chosen
	def exit_if_no_method_selected(self ):
		args=self.parser.parse_args()
		if not self.method:
			print "choose a method for some action... \nchoices:\n    "+"\n    ".join(self.method_choices)
			exit()

	#call this in setup_run to make sure that user doesn't try to modify the setup
	def check_absence_method_options(self):
		if not self.method_arg_group: return
		for i in self.method_arg_group._group_actions:
			for str in i.option_strings:
				if str in self.full_cmdline:
					raise library.RunException("illegal option: %s!\nOptions that change the setup are not allowed in this application. Use setup_target.py to modify the setup first"%str)

	#set refreshed instance of method
	def fresh_method( self ):
		if self.method:
			self.method=self.method.fresh_instance()
		return self.method

	#generate a targetdir from cmdline options and target_name, target_name can be None
	def get_target(self, target_name, rundir=None ):
		args=self.parser.parse_args()
		if target_name:
			target=TargetDir( target=target_name, prefix=args.target_prefix, dir=args.target_dir, rundir=rundir );
		else: #use default target-name 't000_'
			target=TargetDir( dir=args.target_dir, prefix=args.target_prefix, rundir=rundir )
		return target

	#load the options of the setup
	def load_setup(self, setup, verbose=True ):
    #obtain instance of Setup
		if self.hidden_parser:
			method_args=self.hidden_parser.parse_known_args()[0]
			args=self.parser.parse_args()
			self.method.set_args(args)
		else:
			method_args=self.parser.parse_args()
			args=method_args
	  #update options
		self.method.extract_method_options( method_args, self.method_arg_group._group_actions )
		if setup.exists():
			args=self.method.load_options( setup, library.flatten_list( args.remove ), parser=self.parser )
			if verbose: self.method.show("LOADED: Method options from existing setup '%s'..."%setup.name)
		return args

	#modify the options in the setup using the cmdline and store the result
	def modify_setup(self, target, setup ):
		args=self.parser.parse_args()
		self.method.extract_method_options( args, self.method_arg_group._group_actions )
		self.method.upload_files( target )
		target.write_file_database()
		self.method.store_options( setup )
		self.method.show("STORED: method options as new setup '%s'..."%setup.name)

