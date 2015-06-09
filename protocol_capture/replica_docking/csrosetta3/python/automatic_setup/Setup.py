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
import library
import os
from os import path
from os.path import exists
import TargetDir
#Setup class is responsible for the combination of method and target
# it creates a subtree with <method>/<label>
# files that are created out of original input files during target processing (RunGeneration)
# are stored in the setup-root, rather than in the target-root.
# e.g., centroid-mapped constraints
# or upl-generated noQF, QFall files.
# the RunGenerator talks to the TargetDir only through Setup
#
# a setup stores an "option" file which is basically the cmd-line options used to generate it.
# these options are used by the RunGenerator to translate this into a runnable directory


tr = library.Tracer( "Setup" )

class Setup:
	def __init__(self, target, method_name, label, regenerate=False ):

		#method name
		self.method_name=method_name
		#method-variant
		self.label=label
		#protein target
		self.target=target

		#name of setup METHOD_LABEL
		self.name=self.method_name+"_"+label
		self.setup_dir=path.normpath(self.method_name+"/"+label).replace('//','/')
		self.option_library_file=target.dir()+"/"+self.setup_dir+'/options.txt'
		self.regenerate_derived_files=regenerate

		#dictionary for options
		self.options={}
#		self.patches={}

#		library.mkdirp( self.target.dir()+"/"+self.setup_dir )
		self._run_dir = None
		self._generated_files = []
		self._generated_abspaths = {}

	def set_run_dir( self, rundir ):
		self._run_dir = rundir

	def __str__(self):
		return "Setup( %s | %s )"%(self.target.name,self.name)

	#target exists if there is the directory and the option-file
	def exists(self):
		if exists(self.target.dir()+"/"+self.setup_dir):
			return exists(self.option_library_file)
		else:
			return False

	#destructor: if the setup is without any information (now options file) remove the directories again.
	def __del__(self):
		if not self.exists():
			try:
				os.removedirs(self.setup_dir)
			except:
				pass

	#the root of the setup relative to the target-root
	def dir(self):
		return self.setup_dir

	#check if file is in target-database
	#look first in setup sub-tree, then in the general target-tree
	def has_file(self, file, subdir=None, setup_only=False ):
		#print file, subdir
		sfile, dir = TargetDir.subdir_name(file, subdir, 'inputs' )
		if sfile in self._generated_files:
			return True, sfile
		if file in self._generated_files:
			return True, file
		if setup_only:
			return False, ""
		has, file = self.target.has_file( file, subdir=subdir, setup_dir=self.setup_dir )
		if has:
			return has and not self.regenerate_derived_files, file
		has, file = self.target.has_file( file, subdir=subdir )
		return has and not self.regenerate_derived_files, file

	#for files created during RunGeneration
	def create_file(self, file, subdir=None ):
		tr.out("asked to create new file %s"%file)
#		return self.target.add_file( file, subdir=subdir, setup_dir=self.setup_dir, new=True )
		new_file, subdir_str = TargetDir.subdir_name( file, subdir, 'inputs' )
		library.mkdirp( self._run_dir+'/'+subdir_str )
		self._generated_files.append(new_file)
		self._generated_abspaths[new_file]=path.normpath(self._run_dir+'/'+new_file).replace('//','/')
		tr.out("created new_file %s at abspath %s"%(new_file,self._generated_abspaths[new_file]))
		return new_file

	#return CM_TARGETPATH -type path.
	#look first in setup sub-tree, then in the general target-tree
	def cm_path(self, file, subdir=None ):
		tr.out("cm_path: asked for file %s in %s"%(file,subdir))
#		print self._generated_files
		has_self, sfile = self.has_file( file, subdir, setup_only=True )
		if has_self:
			tr.out("have found file as %s..."%sfile)
		else: tr.out("have not found file... relay to TargetDir")
		if has_self:
			return "../"+sfile
		try:
			return self.target.cm_path( file, subdir=subdir, setup_dir=self.setup_dir )
		except:
			return self.target.cm_path( file, subdir=subdir )

	#return absolute paths
	#look first in setup sub-tree, then in the general target-tree
	def abspath(self, file, subdir=None):
		#print file, subdir
		has_self, sfile = self.has_file( file, subdir, setup_only=True )
		if has_self:
			return self._generated_abspaths[sfile]
		elif file in self._generated_files:
			return self._generated_abspaths[file]
		try:
			return self.target.abspath(file,subdir=subdir, setup_dir=self.setup_dir )
		except:
			return self.target.abspath(file,subdir=subdir)

	#print file list (which is handled by target)
	def print_file_list(self):
		print self.target

	#load options from file generates list of triples (cmdline-arg, value(s), is_list(bool))
	#this is called by BasicMethod.load_options
	def load_options( self ):
		all=[]
		file=open(self.option_library_file,'r')
		for line in file:
			if len(line)<1: continue
			if line[0]=='#': continue
			line=line[:-1]
			tags=line.split('|')
			opt=tags[0]
			val=tags[1]
			is_list=False
			if val[0]=='[':
				val=[i.strip().strip("'") for i in val[1:-1].split(',')]
				is_list=True
			all.append((opt,val,is_list))
		all.sort(key=lambda x: x[0])
		return all

	#store options in human-readable format
	def store_options( self, target_options ):
		if len(target_options.keys()):
			library.mkdirp( self.target.dir()+"/"+self.setup_dir )
			file=open(self.option_library_file,'w');
			for opt,val in target_options.iteritems():
				if val!=None:
					file.write("%s|%s\n"%(opt,str(val)))
