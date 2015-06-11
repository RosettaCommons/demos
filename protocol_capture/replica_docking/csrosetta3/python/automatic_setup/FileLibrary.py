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

#FileLibrary provides files (such as flag-files and patches) that are organized like dictionaries:
#e.g., the function of a line is specified by its first word and each key-word appears only once.
#the sequence of lines does not matter.
#this is true for flag-files and patches of rosetta
#key-files are manipulated by provide_file, override, add and remove
#
#additionally the FileLibrary also handles the entry "executable and tagged_Strings
# these should maybe be handled by a different class entirely as they are independent from the files part.
#
#data structure:
# a line_dictionary is ( dic, lines) where dic maps flag->line-number and the respective line can be accessed:
# lines[ dic[flag] ]
#
class FileLibrary:
	def __init__( self, method ):
		self.files = {}
		self.strings = {}
		self.tagged_strings = {}
		self.files["patches"] = {}
		self.files["flags"] = {}
		self.files["others"] = {}
		self.executable = "minirosetta"
		self.tagged_strings["commandline"] = []
		self.method = method # for exceptions

	#make this FileLibrary provide a file <name> (e.g., flags_denovo) under <entry> (e.g., patches, flags... )
	#found on the disk at <path>+<name> or at <source>
	def provide_file( self, entry, path, name, source=None ):
		dic=self.files[ entry ]
		if name in dic:
			raise MethodException(self.method,"cannot provide "+name+" to "+entry+" because it already exists")
		if not source:
			source=path+name
		if not exists( source ):
			raise MethodException(self.method,"cannot provide "+name+" to "+entry+". File not found !!! ")
		dic[ name ] = self._file2dic( source )

	#override individual flags: apply patch to flag-file <name> in category <entry>
	#example for patch:  "-iterative:max_nstruct  -1 -1 0 0 0 0"
	#first word in patch will be matched against existing flags in <name> and these lines will be replaced by <patch>
	def override( self, entry, name, patch ):
		dic=self.files[ entry ]
		if not name in dic:
			raise MethodException(self.method, "cannot override something in "+name+" because it is not present yet");
		line_dic=dic[name]
		tags=patch.split()
		if not tags[ 0 ] in line_dic[0]:
			raise MethodException( self.method, tags[ 0 ]+" not found in "+name+"; cannot override ");
		line_dic[1][ line_dic[0][ tags[0] ] ]=patch

	def get_line( self, entry, name, patch ):
		dic=self.files[ entry ]
		if not name in dic:
			raise MethodException(self.method, "cannot override something in "+name+" because it is not present yet");
		line_dic=dic[name]
		tags=patch.split()
		if not tags[ 0 ] in line_dic[0]:
			raise MethodException( self.method, tags[ 0 ]+" not found in "+name+"; cannot override ");
		return line_dic[1][ line_dic[0][ tags[0] ] ]

	#add to flag_file <name> in <entry>
	#e.g., -rdc med_new.rdc
	#if -rdc not present in file, it just adds a new line "-rdc med_new.rdc"
	#if -rdc is present in file as "-rdc med_old.rdc", it will change the flag to "-rdc med_old.rdc med_new.rdc"
	def add( self, entry, name, patch ):
		dic=self.files[ entry ]
		if not name in dic:
			raise MethodException(self.method, "cannot add to "+name+" because it is not present yet");
		line_dic=dic[name]
		tags=patch.split()
		if not tags[ 0 ] in line_dic[0]: #flag is not yet contained in file
			line_dic[ 0 ][ tags[0] ]= len( line_dic[ 1 ] ) #the number of the last-line +1
			line_dic[ 1 ].append( patch ) #adding the line at the end
		else: #flag has been present, figure out line-number and append to string of that line
			line_num = line_dic[ 0 ][ tags[0] ]
			line_dic[ 1 ][ line_num ] = line_dic[ 1 ][ line_num  ] +" "+" ".join( tags[1:] )

	#remove line from flag-file
	#e.g., remove line='-rdc' will remove line with flag -rdc
	def remove( self, entry, name, patch ):
		dic=self.files[ entry ]
		if not name in dic:
			raise MethodException(self.method, "cannot remove something in "+name+" because it is not present");
		line_dic=dic[name]
		tags=patch.split()
		if not tags[ 0 ] in line_dic[ 0 ]:
			raise MethodException( self.method, tags[ 0 ]+" not found in "+name+"; cannot remove ");
		line_num = line_dic[ 0 ][ tags[0] ]
		del (line_dic[ 1 ])[ line_num ]
		del (line_dic[ 0 ])[ tags[0] ]

	#called by Substitutor to obtain the line-array of a file in the Filelibrary.
	#the caller doesn't have to specify the category
	def read_lib_file( self, filename ):
		for file_dic in self.files.itervalues():
			if filename in file_dic:
				return self._dic2lines( file_dic[ filename ] )
		raise MethodException( self.method, "cannot find "+filename+" in file-library" )

	#add a single tag to a tagged string
	def add_tag_to_string( self, entry, string ):
		self.tagged_strings[ entry ].append( string )

	#in <entry> (e.g., commandline) add more tags (e.g., words) as they are found in string
	#add_string( 'commandline', '-flag1 -flag2' ) will add the tags '-flag1', '-flag2'
	def add_string( self, entry, string ):
		tags=string.split()
		dic=self.tagged_strings[ entry ]
		self.tagged_strings[ entry ]= dic +  tags

	#assembly the string <entry> from its list of tags
	def assemble_string( self, entry ):
		dic=self.tagged_strings[ entry ]
#		print "command-line as tags: ", dic
		str = ""
		for s in dic:
			str = str + " " + s
		return str


	def show( self ):
		print " ========== FILE-LIBRARY ================"
		for (type, file_dic) in self.files.iteritems():
			print type+": { "
			indent1="     ";
			for (file, lines ) in file_dic.iteritems():
				print "%s %s : ["%(indent1,file )
				indent2=indent1+"     ";
				for line in lines[1]:
					print "%s %s"%(indent2, line)
				print "%s ]"%indent1
			print "}"
		for ( type, line ) in self.tagged_strings.iteritems():
			print type+": { "+" ".join(line)+" }"
		print "executable: ", self.executable
		print " ========================================"


	#private methods
	#
	def _file2dic( self, file ):
		#dic[0] : dictionary flag -> line-number
		#dic[1] : lines
		lines = open( file, 'r' ).read().split('\n')
		ll = []
		for l in lines:
			ls=l.split('###')
			if '###' in l and len(ls[0])<1: continue
			ll.append(ls[0])
		lines=ll

		dic = ({}, lines )
		for num,line in enumerate(dic[1]):
			tags=line.split()
			if len(tags) > 0:
#				print tags
				if not tags[0][0]=='#':
					dic[0][ tags[0] ]=num
#		print dic
		return dic

	def _dic2lines( self, dic ):
		return dic[1]


