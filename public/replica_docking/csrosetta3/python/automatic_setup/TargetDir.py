#!/usr/bin/python
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

import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen,mkdir,makedirs
from os import dup2
from operator import add
from math import sqrt
from os import path
from library import MissingInput
import os

import argparse
import sys
import shutil
### toolbox library
import library
from library import StubbedOut
import socket

#class is responsible for a single target
# * it keeps track of the files added to all setups with this target in the list files
# * it also keeps a dictionary with the original paths for files
# * when instantiated a dictionary with absolute paths is generated (i.e., target-libs can be moved on the file-system without depraciating the data)
#
def subdir_name(file,subdir,setup_dir):
	if not subdir and '/' in file:
		subdir=path.dirname(file)
	subdir_str=""
	if setup_dir:
		subdir_str=subdir_str+"/"+setup_dir
	if subdir:
		subdir_str=subdir_str+"/"+subdir
	subdir_str=path.normpath(subdir_str[1:]).replace('//','/')
	new_file=path.normpath(subdir_str+"/"+path.basename(file)).replace('//','/')
	return new_file, subdir_str


class TargetDir:
	def __init__(self, dir=None, target='t000_', prefix=None, overwrite=False, rundir=None ):

#the name of the target
		self.name=target

#overwrite files if they are re-added despite different checksum ?
		self.overwrite=overwrite

		self.files=[]
#figure out the target-directory
		self.fn_target_dir=""
		if dir: #it is a one-off directory
			self.fn_target_dir=dir
		elif prefix: #it is in the targetlibrary
			self.fn_target_dir=path.normpath(prefix+"/"+target).replace('//','/')
		else:
			raise MissingInput("require either directory or directory-prefix to figure out target directory" )

	#	if not path.exists( self.fn_target_dir ):
	#		library.mkdirp( self.fn_target_dir )
		#this prefix is the relative position from the final run-directory if relative paths are used
		self.target_dir_prefix_for_run=""
		self.fn_target_dir=path.normpath(self.fn_target_dir).replace('//','/')
#		print "TARGET_DIR: ", rundir, self.fn_target_dir
		if rundir and not path.isabs(self.fn_target_dir):
			if path.isabs(rundir):
				self.fn_target_dir=path.abspath(self.fn_target_dir)
			else:
				self.rel_path_target_dir=path.normpath(path.relpath(self.fn_target_dir, rundir)).replace('//','/')
				com=path.commonprefix([rundir, self.fn_target_dir ])
				if com: print "Common-Prefix: ", com
				target_dir_alone=self.fn_target_dir.replace(com,"")

				self.target_dir_prefix_for_run=self.rel_path_target_dir.replace(target_dir_alone,"")
				print "relative path from run-dir to input dir", self.rel_path_target_dir, self.target_dir_prefix_for_run
#				com=path.commonprefix([rundir, self.fn_target_dir ])
#				If com: print "Common-Prefix: ", com
#				target_dir_alone=self.fn_target_dir.replace(com,"")
#				rundir_alone=rundir.replace(com,"")
#				print "correcte paths", target_dir_alone, rundir_alone


		#keep track of the original_paths  -- keep the hostname
		self.hostname=socket.gethostname()
		self.original_paths={}

		#secondary data -- generate from stored file-list at each instantiaten
		self.absolute_paths={}

    #the master list of files stored in this directory (with add_files)
		self.directory_file=self.fn_target_dir+'/file_directory.txt'
		if path.exists( self.directory_file ):
			self.read_file_database()
#			print "found target directory at %s ..."%(self.fn_target_dir )

			#root of target directory
	def dir(self):
		return self.fn_target_dir

	#check if there are files stored in this target -- this might be wrong... what if there is setup that soley uses options...
	# but it is hard to see working setups that don't require any input files
	def exists(self):
		return len(self.files)

	#add files to the targets file library, use new=True to create new files
	def add_file(self, file, new=False, subdir=None, setup_dir=None ):

		# if this is not a new file, we require that there is a physical file we can copy
		if not new and not path.exists( file ):
			raise MissingInput( "Cannot find file %s"%file)

		# if this is the first file we add, make a directory
		if not path.exists( self.fn_target_dir ):
			library.mkdirp( self.fn_target_dir )

		#figure out relative filename
		new_file, subdir_str = subdir_name( file, subdir, setup_dir )
#		print "TARGET: ", new_file, subdir_str, subdir, setup_dir

		#generate subdirs if necessary
		library.mkdirp( self.fn_target_dir+"/"+subdir_str )

		#check if file is already known
		if self.overwrite: overwrite_str=' going to overwrite'
		else: overwrite_str =''
		if new_file in self.files:
			if not new:
#				print "check ", new_file, file
				if library.diff_files_md5( self.abspath(new_file), file):
					print "%-50s was present in previous setup, md5 checksum OK ..."%(new_file)
				else:
					print "%-50s was present in previous setup, md5 checksum FAIL, %s ..."%(new_file,overwrite_str)
					if not self.overwrite:
						raise library.InconsistentInput("%s was present in previous setup, md5 checksum FAILED...\n"%new_file+
																					"check your input...\n"
																					"use -overwrite to overwrite files\n"+
																					"CAREFUL: overwriting of files affects all Setups of Target %s (all labels, all methods)\n"%self.name+
																					"better to use different filename unless this is not a correction of an errorneous file")
		else:
			#add this file as it hasn't been present already
			self.files.append(new_file)
			self.absolute_paths[new_file]=path.normpath(self.fn_target_dir+"/"+new_file).replace('//','/')
			if file[0]=='/':
				self.original_paths[new_file]=self.hostname+":"+file
			else:
				self.original_paths[new_file]=self.hostname+":"+path.normpath(os.getcwd()+"/"+file).replace('//','/')

			#copy the data, if it hasn't been a new file creation 'new=True'
		if not new:
			shutil.copy(file, self.absolute_paths[new_file]);

		#return the filename relative to the root of TargetDir
		return new_file

	#return the root of the TargetDir
	def dir(self):
		return self.fn_target_dir

	#figure out the filename and subdir within TargetDir
	# say we look for nmr/cs.tab
	# subdir_name(cs.tab, nmr, None ) -> nmr/cs.tab, nmr
	# subdir_name(nmr/cs.tab, None, None ) -> nmr/cs.tab, nmr
	# subdir_name(nmr/cs.tab, None, abrelax/standard ) --> abrelax/standard/nmr/cs.tab, abrelax/standard/nmr
	# subdir_name(cs.tab, nmr, abrelax/standard ) -->  abrelax/standard/nmr/cs.tab, abrelax/standard/nmr
	# subdir_name(abrelax/standard/nmr/cs.tab, None, None ) --> abrelax/standard/nmr/cs.tab, abrelax/standard/nmr

	#absolute path of file
	def abspath( self, file, subdir=None, setup_dir=None ):
		try:
			sfile, dir = subdir_name( file, subdir, setup_dir )
			if sfile in self.files:
				return self.absolute_paths[sfile]
			else:
				return self.absolute_paths[file]
		except KeyError:
			self.raise_missing_file_error( file )

		#check if file is present in TargetDir
		#ATTENTION: returns also the filename as which it was found
		# return bool, str
		#logic: either file is already the fully specified name relative to the root of TargetDir
		# or file is build up using subdir and setup_dir.
		# has,file = has_file( 'cs.tab', subdir=nmr ) -> true, nmr/cs.tab
		# has,file = has_file( 'nmr/cs.tab', subdir=nmr ) -> true, nmr/cs.tab
	def has_file( self, file, subdir=None, setup_dir=None, exception=False ):
		sfile, dir = subdir_name(file, subdir,setup_dir )
#		print file,sfile, self.files, '\n\n'
		if exception:
			if sfile in self.files:
				return sfile
			if file in self.files:
				return file
			self.raise_missing_file_error( file, sfile )
		else:
			if sfile in self.files:
				return True, sfile
			if file in self.files:
				return True, file
			return False, ""

	#return path with CM_TARGETPATH, logic similar to hasfile
	def cm_path( self, file, subdir=None, setup_dir=None ):
#		sfile, subdir=subdir_name(file, subdir, setup_dir)
		file = self.has_file( file, subdir, setup_dir, exception=True );
		return path.normpath("$CM_TARGETPATH/"+file).replace('//','/')

	def cm_dir( self, levels ):
		prefix=""
		if self.target_dir_prefix_for_run:
			prefix="../"*levels+self.target_dir_prefix_for_run
		print "CM_DIR: ", prefix+'/'+self.dir()
		return path.normpath(prefix+'/'+self.dir()).replace('//','/')

	#printable summary of the Target including a list of the files stored in the database
	def __str__( self ):
		format="*"*60+'\n'
		format+="Target %20s   ---  Path: %s\n"%(self.name,self.fn_target_dir)
		format+='-'*60+'\n'
		for f in self.files:
			format=format+"->%s<-\n"%f
		format+='*'*60
		return format

	#printable just the name of Target
	def __repr__( self ):
		return self.name

	#destructor,
	#  1) write file database
	#  2) remove directory if empty
	def __del__(self):
		if len(self.files):
			self.write_file_database()
		else:
			try:
				os.removedirs(self.dir())
			except:
				pass

	#write files and original_paths
	def write_file_database(self):
		if path.exists(self.directory_file ):
			shutil.copy(self.directory_file, self.directory_file+".backup" )

		dir=open(self.directory_file,'w')
		for file in self.files:
			dir.write( "%-60s %s\n"%(file, self.original_paths[file] ))

	#read files and original_paths
	#generate absolute_paths with current root of TargetDir
	def read_file_database(self):
		if not path.exists(self.directory_file):
			return
		dir=open(self.directory_file,'r')

		self.files=[]
		self.absolute_paths={}
		self.original_paths={}

		for line in dir:
			if len(line) < 1: continue
			if line[0]=='#': continue
			col=line.split()
			self.files.append(col[0])
			self.original_paths[col[0]]=col[1]
			self.absolute_paths[col[0]]=(self.fn_target_dir+"/"+col[0]).replace('//','/')

	#report error
	def raise_missing_file_error( self,file, alt=None ):
		if alt:
			raise library.ProgramError("cannot find file %s (also looked for %s) in Target %s\n(usually added by <method>/options.py during setup)"
																 %(file,alt,self.name))
		else:
			raise library.ProgramError("cannot find file %s in Target %s\n(usually added by <method>/options.py during setup)"%(file,self.name))
