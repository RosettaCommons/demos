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
from os.path import basename,dirname
from os.path import exists
from os.path import splitext

from sets import Set
import sys
import shutil
### toolbox library
from library import *

### bla=@file  bla=@$CM_XX $( bla -x $CM_XX) `$CM_EXEC -x ` bla=@@file

## strategy;
# tokenize with all whitespaces but '\' and run the SUBS --> replaces $CM_XXX but not \$CM_XXX
# tokenize with whitespaces but not '@','.','!', and '/' which would keep
# filenames and their @@ in place.
# substitute filenames that are relative and start with fn_rundir/ with ../
# start substitution on new files
# tokenizing is done recursively one special character at a time
# when all characters are used the "action" can be variable:
# i.e., substitute_CM_tag or substitute_file

tr = Tracer( "Substitutor" )

class Substitutor:
	def __init__( self ):
		pass

	def whitespace( self ):
		pass

	def apply( self, str ):
		pass

class File:
	def __init__(self, file, fn_rundir):
		self.input_fn=file #where it is now
		self._output_fn=file  #where to write this while setup is running
		self.tag_in_file="" #where to find this when rosetta is running
		if not "/" in file: #assume to find this in flag_library
			pass#			self._output_fn=fn_rundir+"/"+self.input_fn
		else:
			#self._output_fn=fn_rundir+"/"+basename(self.input_fn)
			if not exists( self.input_fn ):
				if self.input_fn[0:2]==".." and exists (self.input_fn.replace("..",fn_rundir,1)):
					print 'found here'
					self.input_fn=self.input_fn.replace("..",fn_rundir,1)
				else:
					raise MissingInput( "file %s not found"%self.input_fn)
		self.tag_in_file=basename(self.input_fn)

	def output_fn( self, destdir ):
		if not "/" in self._output_fn: #assume to find this in flag_library
			return destdir+"/"+self._output_fn
		else:
			return destdir+"/"+basename(self._output_fn)

	def __str__(self):
		return self.input_fn

	def __eq__(self, other):
		return self.input_fn == other.input_fn

	def __hash__(self):
		return hash('%s'%self)

	def read_lines(self, method):
		lines=[]
		if not "/" in self.input_fn:
			lines=method.read_lib_file( self.input_fn )
			tr.out("opened %s via method-lib"%self.input_fn )
			tr.out("\n".join( lines ) )
#			print "open via method-lib"
		else:
			try:
				lines=open(self.input_fn,'r').read().split('\n')
				tr.out("opened %s directly from disc"%self.input_fn )
				tr.out("\n".join( lines ))
#				print "found file in target-lib: %s"%self.input_fn
			except IOError:
				raise MethodException(method, "file %s not found"%self.input_fn)
		return lines

class SubstituteFiles( Substitutor ):

	def __init__( self, rundir ):
		self.whitespace=[' ','"','`',"'",'=']
		self.file_list=[]
		self.rundir=rundir

	def whitespace( self ):
		return self.whitespace

	def apply( self, tag ):
		if len( tag )<2:
			return tag
		AT='@'
		path_corrected = self.fix_relative_path( tag )
		if tag[0] != '@':
			return path_corrected
		tag=tag[1:]
		if tag[0]=='@':
			AT=''
			tag=tag[1:]
					#now we xxxxxxhave @ or @@ at beginning of tag --> file
		new_file = File( tag, self.rundir )
		self.file_list.append( new_file )
		return AT+new_file.tag_in_file

	def fix_relative_path( self, tag ):
		#	print tag+"-->"+tag

		#this is obsolete:
		#if "/" in tag and not path.isabs(tag):
		#	tag=tag.replace(self.rundir,"..")
		return tag


 	def substitute_file(self, method, filename):
		lines=[]
		if not "/" in filename:
			lines=method.read_lib_file( filename )
			output_fn=self.fn_rundir+"/"+filename
#			input_fn=method.path+"/"+filename
#			print "open via method-lib"
		else:
			input_fn=filename
			filename=basename(filename)
			output_fn=self.fn_rundir+"/"+filename
#			print "open as file"
			try:
				if exists( input_fn ):
					lines=open(input_fn,'r').read().split('\n')
				elif input_fn[0:2]==".." and exists (input_fn.replace("..",self.fn_rundir,1)):
					#	try to replace leading .. with fn_rundir
					alternative_input_fn=input_fn.replace("..",self.fn_rundir,1)
#					print "found file in target-lib: %s"%alternative_input_fn
					lines=open(alternative_input_fn,'r').readlines()
		#				print lines
			except IOError:
				raise MethodException(method, "file %s not found"%input_fn)


#		print lines
#		print "substitute flag file %s --> %s"%(filename, output_fn)
		out_file=open(output_fn,'w')

		for line in lines:
			out_file.write( self.substitute(method, line) )
		rel_path=""
		if not self.all_flags_in_rundir:
			rel_path="../"
		return rel_path+basename(output_fn)


class SubstituteCMTags( Substitutor ):
	def __init__( self, SUBS ):
		self.whitespace=[' ','/','"','.','`',"'",'=','!','@']
		self.SUBS=SUBS

	def apply( self, tag ):
		if len( tag )<5:
			return tag
		if tag[0:4]=="$CM_":
			tag=tag[1:]
		else:
			return tag

		subst_tag = self.SUBS[tag]
		return subst_tag

class Tokenizer:
    def __init__( self, substitutor ):
 		  self.substitutor=substitutor
		  self.whitespace=substitutor.whitespace

	 def subst_token( self, str ):
		 return self.substitutor.apply( str )

	 def has_any_whitespace( self, str ):
		 for i,c in enumerate(self.whitespace):
			 if c in str:
				 return i
		 return None

    def apply( self, token_in, whitespace_index=0 ):

		 #done all whitespace ? run the substitutor on the final token.
		 if whitespace_index >= len( self.whitespace ):
			 subst_tag = self.subst_token( token_in )
			 ci = self.has_any_whitespace( subst_tag )
			 if ci:
				 subst_tag=self.apply( subst_tag, ci )
	#		 tr.out("return string from tokenizer: %s"%subst_tag )
			 return subst_tag

		  #still more whitespaces to take care of
		 c = self.whitespace[ whitespace_index ]
		 if c in token_in:
			 token_out =""
	#		 tr.out("tokenize %s with whitespace: %c"%(token_in, c) )
			 tags=token_in.split( c )
			 for t in tags:
				 token_out = token_out + self.apply( t, whitespace_index+1 ) + c
			 return token_out[:-1] #remove superfluous c at end
		 else:  #c is not present, skip this level of tokenizing
			 return self.apply( token_in, whitespace_index + 1 )



