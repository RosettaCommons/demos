#!/usr/bin/env python2.7
##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
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
from glob import glob
#from sys import argv,stderr,exit
#import sys
import subprocess
from os import dup2,path
from operator import add
from math import sqrt
import os
from os.path import basename,dirname
from os.path import exists
from os.path import splitext
from os import system
from Job import Job
from sets import Set
from Substitutor import *
import sys
import shutil
### toolbox library
from library import *
from Setup import Setup

tr = Tracer( "RunGenerator" )

#RunGenerator is responsible to synthesize the information of the method and target into a runnable set of files
#it uses the Substituor to replace placeholders like "CM_EXECUTABLE" with actual files and correct pathnames.
# it uses "BasicMethod" to
class RunGenerator:

	def __init__(self, rundir, flag_lib, platforms, bindir, database=None, extension=".default.linuxgccrelease", subdir_level=1 ):
		if rundir[-1]=="/": rundir=rundir[0:-1]
		self._fn_rundir = rundir
		self.fn_flag_lib = flag_lib
		self.platforms=platforms
		self.all_flags_in_rundir=True #have flags in actual run-dir or in top-level directory : @../flag or @flag
		self.file_list=Set()
		self.rerun=False
		if not database:
			raise library.MissingInput('Need database: either specify ROSETTA3_DB environment variable or use -database option')
		if database[0]!="/":
			db_prefix="../../"
		else:
			db_prefix=""
		self.database=db_prefix+database
		self.bindir=bindir
		self.extension=extension
		self._extra_test=True
		self.subdir_level=subdir_level
		self._substituted_files=Set()

	def rundir(self):
		return self._fn_rundir
		#generate the initial directory structure: rundir/run and rundir/test
	def generate_dirs(self, overwrite=False):
		if exists(self._fn_rundir):
			if not overwrite:
				raise RunException("directory exists already: %s ..."%self._fn_rundir)
			else:
				system('rm -Rf %s'%self._fn_rundir)
		#mkdirp( self._fn_rundir );
		mkdirp( self._fn_rundir+"/run" )
		if self._extra_test: mkdirp( self._fn_rundir+"/test")

		#figure out the name and path of the executable and store in SUBS dictionary
	def setup_executable(self, method):
		exe=method.executable()
		full_path_exe=self.bindir+"/"+exe
		if full_path_exe[0]!="/":
			full_path_exe="../../"+full_path_exe
		self.SUBS["CM_EXECUTEABLE"] =  full_path_exe#make this another option
		self.SUBS["CM_EXEC_EXT"] = self.extension

		#this is a flag-file generated directly from the options instead of going via flag-library
		# maybe obsolete ?
		# things like -nstruct or -in:file:rdc goes here.
		# this is the place for the most input files, if they are not introduced via flag-library and CM_XXX substitutions (like frags)
		# this method calls 'make_target_flags', the method specific code to setup flags
	def setup_target_flags(self, method, setup ):
		self.targetflags=setup.create_file("flags_"+method.name)
		print self.targetflags
		str_fn=setup.abspath(self.targetflags)
		print "RunGenerator: generate targetflags %s as %s"%(self.targetflags, str_fn)

		flag_file=open(str_fn,'w')
		flag_file.write("#this flag-file is automatically generated\n\n") #TEMPORARY
		method.make_target_flags(self, setup, self.targetflags, flag_file, self.SUBS )
		self.SUBS["CM_FLAGFILE"] = str_fn   #absolute paths for flagfiles, so they get substituted right

		#setup the full commandline (everything behind the executable)
	def setup_commandline( self, method ):
		cmdline=method.commandline()
		print "generic commandline: ",cmdline
#		cmdline=self.substitute(method, cmdline)
#		print "CMDLINE: ",cmdline
		self.SUBS["CM_COMMANDLINE"] = cmdline

	def copy_files( self, target_dir ):
		print 'copy flag-files to %s'%target_dir
		for file in self.file_list:
#			print "copy file %s to %s"%( file, target_dir )
			shutil.copy( self._fn_rundir+"/run/"+file, target_dir )

	def copy_final_files_to_rundir(self):
		flag_list=open(self._fn_rundir+"/run/flag_list.txt",'w')
		self.file_list.add('flag_list.txt')
		for file in self.file_list:
			flag_list.write("%s\n"%(file))
		flag_list.close()

		if self._extra_test:
			for dir in ['test']:
			   self.copy_files( self._fn_rundir+"/"+dir )

		#create a script for copying in case somebody creates run_2 or so
		script=open(self._fn_rundir+"/run/copy_flag_files.sh",'w')
		os.fchmod(script.fileno(),0o755)
		script.write("mkdir -p $1\n")
		for file in self.file_list:
			script.write("cp %s $1\n"%(file))
		script.write("cp *.job $1\n")
#		script.write("cp flag_list.txt $1\n")
		script.write("cp copy_flag_files.sh $1\n")

		setup_cmd=open(self._fn_rundir+"/setup_command.csh",'w')
		args=sys.argv;
		dels=[]
		if '-overwrite' in args:
			del args[args.index('-overwrite')]
		if '-dir' in args:
			i=args.index('-dir')
			dels.append(args[i])
			dels.append(args[i+1])
			del args[i+1]
			del args[i]
		setup_cmd.write(" ".join( args )+"\n")
		setup_cmd.write("# "+" ".join( dels )+"\n");

	def generate_flags(self, setup, method ):
		if not method:
			raise MissingInput("requires definition of a method");

		self.flag_list=Set()
		assert issubclass( setup.__class__, Setup )
		setup.set_run_dir( self._fn_rundir )
		print "generate run from ", setup
		self._generate_subs( setup )
		self.setup_executable( method )
		self.setup_target_flags( method, setup )
		self.setup_commandline( method ) #this will generate flags via the substitutions
		tr.out("SUBS: \n"+str(self.SUBS))
		self.generate_jobscripts( method )
		self.copy_final_files_to_rundir()
#		self.generate_template_dir_and_store_cmd( setup )

#### Substitution Syntax:
##  "@@file" --> substitute zu "file" add "file" to list of config files
##  "@file" --> substitute zu "@file" add "file" to list of config files
##  "$CM_XX" --> substitute $CM_XX according to substitution list e.g., CM_NATIVE, CM_METHODPATH, CM_FASTA
##  "@$CM_XXX --> subsitute $CM_XX according to substitution list to YYY  and add "YYY" to list of config files
##  "###" --> line will not be transferred
##  "!XXX --> substitute zu "XXX" and continue with rule 1 but do not add files to list of config files if they come up
##     I fogot what the !XXX is good for... doesn't seem to be used currently
##  "\@" escape: treat @ as normal character

	def generate_jobscripts( self, method ):
		for platform in self.platforms:
#			print "generate job-script: %s"%platform
			filename=basename( platform )
			output_fn=self._fn_rundir+"/run/"+filename
#			print "open as file"
			if "test." in filename:
				output_fn=self._fn_rundir+"/test/"+filename
				if not self._extra_test: continue
			try:
				tr.out("substitute job-script %s --> %s"%(filename, output_fn))
				in_file=open( platform,'r')
				out_file=open(output_fn,'w')
				out_file.write( self.substitute( method, in_file.read() ))
			except IOError:
				raise MethodException(method, "file %s not found"%input_fn)

	def substitute(self, method, str):
		subst_lines=""
		lines=str.split('\n')

		for line in lines:
			l=line.split('###')
			if '###' in line and len(l[0])<1:
				continue
			line=l[0]

			cm_tokens=Tokenizer( SubstituteCMTags( self.SUBS ) )
			file_substitutor = SubstituteFiles( self._fn_rundir )
			file_tokens=Tokenizer( file_substitutor )

			subst_line  = cm_tokens.apply( cm_tokens.apply( cm_tokens.apply( line ) ) ) #maximum of level 3 CM -> CM -> CM -> real
			tr.out("after substituting tokens: "+subst_line)
			subst_line  = file_tokens.apply( subst_line )
			tr.out("after substituting files:  "+subst_line )
			if len(file_substitutor.file_list): tr.out("files found: -->"+"<-->".join( [ file.input_fn for file in file_substitutor.file_list  ])  +"<--")
			subst_lines = subst_lines+subst_line+"\n" #eat whitespace at end and add the \n
		#	if len(file_substitutor.file_list): print ", ".join(["%s"%x.input_fn for x in file_substitutor.file_list])
			for file in file_substitutor.file_list:
				if file in self._substituted_files: continue
				self._substituted_files.add(file)
				self.file_list.add( file.tag_in_file )
				tr.out("substitute flag file %s --> %s"%(file.input_fn, file.output_fn( self._fn_rundir+'/run' ) ) )
				print "substitute flag file %s --> %s"%(file.input_fn, file.output_fn( self._fn_rundir+'/run' ) )
				lines=file.read_lines( method )
				out_file=open(file.output_fn( self._fn_rundir+'/run' ), 'w')
				for line in lines:
					out_file.write( self.substitute( method, line) )

		tr.out("return from substitute for string beginning with %s"%str[0:min(100, len(str) )])
		return subst_lines

	def add_subst(self, key, val ):
		self.SUBS[key]=val;

	def _generate_subs(self, setup ):
		#def raise_except(i):
		#	raise RunException(" this method with the chosen options cannot be restarted (that is one cannot use flag -use_existing_setup )")

		self.SUBS = {
#			"CM_FRAGMENTS3" : target_lib.fn_frags3,
#			"CM_FRAGMENTS9" : target_lib.fn_frags9,
#			"CM_FASTA" :      target_lib.fn_fasta,
			"CM_TARGET" :     setup.target.name,
			"CM_TARGETPATH" : setup.target.cm_dir(self.subdir_level), #the 1 because we go one level down from -dir (for the target).
			"CM_METHODPATH" : self.fn_flag_lib,
			"CM_ROSETTA_DATABASE" : self.database,
			"CM_BINPATH"     : self.bindir,
			"CM_SCRIPTPATH"  : path.split( __file__ )[0],
			"CM_AUTO_NSTRUCT": "",
			}
#		if target_lib.fn_native: self.SUBS["CM_NATIVE"]= target_lib.fn_native
		self.SUBS["CM_RUNDIR"]=self._fn_rundir+"/run"

		#if self.rerun:
		#	self.SUBS["CM_NO_RESTART"]=raise_except
		#else:
		#	self.SUBS["CM_NO_RESTART"]=lambda x: ""
