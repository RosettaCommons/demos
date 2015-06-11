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

from os.path import exists
from library import MethodException
import library
from FileLibrary import FileLibrary
import argparse
import sys

#this class is responsible for a method which represents a way to generate rosetta-flags to run a certain protocol.
#inherit from this class to implement your own methods in ../flag_library tree.
# examples: .. see denovo_options, rasrec_options, autoNOE/options.py, abrelax/options.py, rasrec/options.py
#
# Method define their own options and can declare subdirectories to handle certain file options
# e.g., all nmr_data -rdc -cs -peaks etc, can go into nmr_data/
# fasta and fragments goes into fragments/
#
# methods maintain a file-library which is mainly responsible for the flag-files.
# in "setup_file_library" files are provided and patches to the files can be made
#
class BasicMethod:
	def __init__(self,name,path):
		#name of the method -- should coincide with directory name in ../flab_library
		self.name=name
		self.path=path

		#this gets populated by setup_file_library of derived classes
		self.file_library = FileLibrary( self )
		self.file_library_is_initialized = False

		#these get populated by constructors of the derived classes
		self.non_file_options=[]
		self.option2dir={}

		#options where the same file might be given twice  (like -shifts for autoNOE )
		# not super happy about this construction... but no better idea currently
		self.double_file_options=[]

		#these get populated by extract_options and load_options
		self.target_options={}
		self.options={}

		#results from a previous run (e.g., phaseI in autoNOE, or centroid-only in docking and abrelax)
		#can be provided in this directory.
		#the method is responsible to find the correct file (assuming that the run in the respective directory has also been carried out with a
		#compatible method
		self._input_dir = None

		#never acces these: use get_args()
		self._args = argparse.Namespace()
#		print "initialized method %s"%self.name

  	#a squeaky clean instance of the same (derived) type
	def fresh_instance(self):
		return self.__class__(self.name, self.path)

	#a breakdown of the options
	def show(self, description):
		print "\n"+"="*60
		print description
		print "- -"*20
		print self
		print "-"*60+"\n"

	#message to display after setup_run is complete ( can containt instructions to the user... )
	def motd(self, rundir):
		print "\nMethod ",self.name, " has been setup in %s ...\nEnjoy!\n\n"%rundir

	#some details where this Method's python code is found
	def module_description(self):
		return "%s method module from %s"%(self.name, self.path)

	#show a breakdown of target_options
	def __str__(self):
		ss="Method: %s\nChosen Options: \n"%self.name
		for opt,val in sorted(self.target_options.iteritems()):
			if val==None: continue
			if isinstance(val,basestring):
				ss+="   %s %s\n"%( opt, val )
				continue
			try:
				ss+="   %s %s\n"%( opt, " ".join(val))
			except:
				ss+="   %s %s\n"%( opt, val )
		return ss[:-1]

	#this is called when the user calls setup_run with -input
	# rundir is set to the directory where the run is setup (for information only)
	def set_input_dir( self, str, run_dir ):
		self._input_dir = str

	def input_dir( self, str ):
		return self._input_dir

	def substitute(self,line):
		#probably obsolete
		return line

#	def set_target(self,target):
#		assert issubclass( target.__class__, TargetDir )
#		self.target=target

  #virtual function for derived classes to setup the file-library
	#make sure to always also call the parent-class function
	def setup_file_library( self ):
		self.file_library_is_initialized=True

	#virtual function for derived classes
	#make sure to always also call the parent-class function
	def make_target_flags(self, run, setup, filename, flags, subs ):
		pass

	#return the executable that belongs to the method
	def executable( self ):
		assert self.file_library_is_initialized
		return self.file_library.executable

	#return command-line that belongs to this method, as specified by the file-library
	def commandline( self ):
		if not 	self.file_library_is_initialized:
			raise library.ProgramError("FileLibrary is not initialized. Call setup_file_library first ! ")

		return self.file_library.assemble_string( "commandline" )

	#called by Substitutor to read a file fro mthe file-library
	def read_lib_file( self, filename ):
		return self.file_library.read_lib_file( filename )

	#extract the options that belong to the method from the cmdline
	# we now this because method options are grouped into group_actions (a parser)
	def extract_method_options( self, args, group_actions ):
		self.options={}
		self.nargs={}
		for i in group_actions:
			#print i.dest, getattr(args,i.dest)
			val=getattr(args,i.dest)
			val=library.flatten_list(val)
			self.nargs[i.dest]=i.nargs
#			print "flattended: ",val
			if len(val)==1:
				val=val[0]
			self.options[i.dest]=val

	#upload files specified in the cmd-line to the target
	def upload_files( self, target_dir ):
		upload_exceptions=[]
#		assert issubclass( target_dir.__class__, TargetDir )
#		print self.target_options
		#loop over all previously extracted options (extract_method_options)
		for opt,val in self.options.iteritems():
			# do nothing if this is marked as non_file_option
			if opt in self.non_file_options:
				self.target_options[opt]=val
				continue

			#if a potential file option, figure out if a subdir is declared
			subdir=None
			if opt in self.option2dir:
				subdir=self.option2dir[opt]

			if library.obj_is_list( val ):
				for file in val:
					has,new_file=target_dir.has_file( file )
#					print 'BASIC_METHOD: ',file,has,new_file
					if not has:
						try:
							new_file=target_dir.add_file(file,subdir=subdir)
						except library.LibException:
							upload_exceptions.append(sys.exc_info()) #process these later -- use exc_info to keep traceback
					else:
						#print "%-50s transfered from previous setup"%(file)
						pass

					#add file to target_options for use in RunGeneration
					if opt in self.target_options:
						#existing entry -- append to list unless double
						if not new_file in self.target_options[opt] or opt in self.double_file_options:
							self.target_options[opt].append(new_file)
					else: #first entry start new list
						self.target_options[opt]=[new_file]
			else: #opt is not a list-option
				if val != None:
					has,new_file=target_dir.has_file( val )
#					print 'BASIC_METHOD: ',val,has,new_file
					if not has: #new file, not a list-option, make entry in target_optons and add to target
						try:
							self.target_options[opt]=target_dir.add_file(val, subdir=subdir)
						except library.LibException:
							upload_exceptions.append(sys.exc_info()) #process these later -- use exc_info to keep traceback
					else:
						#print "%-50s transfered from previous setup"%val
						self.target_options[opt]=val
						pass
		if len(upload_exceptions):
			#this 3-way call to raise is to keep the trace-back intact... maybe not super necessary here...
			raise upload_exceptions[0][0],upload_exceptions[0][1],upload_exceptions[0][2]

	#store current option-settings to setup
	def store_options( self, setup ):
		setup.store_options( self.target_options )

	#load option-settings from setup, suppress options if specified (-remove)
	#if parser specifie set defaults of cmdline-options with the stored-settings
	#this allows to overwrite stored-settings with cmd=line options
	def load_options( self, setup, suppressed_options=[], parser=None):
		if not self.target_options:
			self.target_options={}

    #nothing to load ?
		if not exists( setup.option_library_file ):
			if parser:
				return parser.parse_args()
			else:
				return

		#load raw data as list of (opt, val, is_list) tuples
		data=setup.load_options()
		for opt,val,is_list in data:
			if not opt in self.options: #valid-option ?
				print "option '%s' from setup %s not recognized -- ignored"%(opt,setup.name)
				continue
			if opt in suppressed_options: #suppressed option ?
				print "option '%s' from setup %s is removed"%(opt,setup.name)
				continue
			if not is_list and self.nargs[opt]=='*': #what with options that have nargs=2 ?
				val = [val]

			#set this option as default
			kwargs={opt:val}
			if parser: parser.set_defaults(**kwargs)

			#add to target_options as list or as single item
			if is_list and opt in self.target_options:
				self.target_options[ opt ] += val #list concatenation
			else:
				self.target_options[ opt ] = val


		if parser:
			#		print 'defaults: ',parser._defaults
			args=parser.parse_args()
			return args
		else: return

	#get namespace with method-options only
	def get_args(self):
		#args=argparse.Namespace()
		args=self._args
		for i in self.options.keys():
			setattr(args,i,None)

		for i,j in self.target_options.iteritems():
			setattr(args,i,j)
		return args

	def set_args(self, args):
		self._args=args;
