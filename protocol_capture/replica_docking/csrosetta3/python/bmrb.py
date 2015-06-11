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
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-

import subprocess
from os import dup2,path
from os.path import exists
from operator import add
from math import sqrt
from os.path import *
import argparse
import sys
import copy
import shutil
import string
### toolbox library
#import library
from os import mkdir , makedirs

from warnings import *

import traceback

#from Bio.PDB.PDBParser import *
#from Bio.PDB import PDBIO
#from numpy import square

from library import mkdirp
import library

class BMRB_Type :
	def __init__( self, tags ):
		self.tags=tags
		self.type=int

class BMRB_Int( BMRB_Type ):
	def __init__( self, tags ):
		BMRB_Type.__init__( self, tags )
		self.type=int

class BMRB_Str( BMRB_Type ):
	def __init__( self, tags ):
		BMRB_Type.__init__( self, tags )
		self.type=str

class BMRB_Float( BMRB_Type ):
	def __init__( self, tags ):
		BMRB_Type.__init__( self, tags )
		self.type=float


_tag_dictionary = {'ID': BMRB_Int(['_Atom_chem_shift.ID','_Atom_shift_assign_ID']),
         'RESID': BMRB_Int(['_Atom_chem_shift.Comp_index_ID','_Residue_seq_code']),
         'RESNAME': BMRB_Str(['_Atom_chem_shift.Comp_ID','_Residue_label']),
         'ATOMNAME': BMRB_Str(['_Atom_name','_Atom_chem_shift.Atom_ID']),
         'SHIFT': BMRB_Float(['_Chem_shift_value','_Atom_chem_shift.Val','_Chem_shift_value']),
         'SIGMA': BMRB_Float(['_Chem_shift_value_error','_Atom_chem_shift.Val_err']) }

cs_loop_cols=['ID','RESID','RESNAME','ATOMNAME','SHIFT','SIGMA']

#after reading the BmrbFile content is organized in frames of name <frame> which are organized in self._frames
#according to their category <frame-category>. A frame goes from 'save_<frame>' to 'save_'
#according to this entry:
#save_<frame>
#   _Saveframe_category      <frame-category>
# .... content
#save_
#From each frame we currently store only the 'loop', that is lines between 'loop_' and 'stop_'. The beginning of a
#loop has entries starting with '_XXX' which give the column names.
#we store each loop of a frame as object of class 'Loop' which has the data-members cols (the column headers) and data
#which is just a list with an entry per line. Each line-entry is a list of tags for this line.
#to extract data from a frame we use 'process_frame( <frame_category>, columns )'
#the columns are translated using the _tag_dictionary which takes care of ambiguous nameing found in different BMRB styles.
#the output is a dictionary {'colname1': [data1, data2, data3, ... , dataN ], 'colname2':[data1, data2, ..., dataN] }
#if multiple frames of the same category have fitting columns these will be appended to the dictionary...

#reads a BmrbFile into _frames
class BmrbFile:
	class Loop:
		def __init__(self, cols, data):
			self.cols=cols
			self.data=data
		def __str__( self ):
			return "[ "+", ".join(self.cols)+" ]"

		def __repr__( self ):
			return self.__str__()

	class Frame:
		def __init__( self, name, fields, loops ):
			self.name=name
			self.fields=fields
			self.loops=loops
		def __repr__( self ):
			str="Frame %s: \n"%self.name
			for k,f in self.fields.iteritems():
				str=str+"%s: %s\n"%(k,f)
			str=str+"and %d loops\n"%len(self.loops)
			return str

	#BmrbFile
	def __init__(self, file):
		self._frames={}
		self.parse_file( open(file,'r') )

	#find next save_<frame>
	#return the name, i.e., <frame> and the category <frame-category>
	def Next_Save_Frame( self, file ):
    line=file.readline();
		name=''
    while line:
			tags=string.split(line)
			#        print tags
			if len(tags)>0 and len(tags[0])>=5 and tags[0][:5]=='save_':
				name=tags[0][5:]
				while len(tags)<1 or tags[0] != '_Saveframe_category':
					line=file.readline();
					tags=string.split(line)
				return tags[1], name
			line=file.readline();
    return 'NO_CATEGORY', 'NO_NAME'

	#store the loops of the current frame as self.Loop
	def capture_loops( self, file ):
		loops=[]
		fields={}
		multi_line_field=None
		line='something to start while loop';
		col_nr=-1;
		while line:
			line=file.readline();
			tags=string.split(line)
			if len(tags)>0 and tags[0]=='loop_':
        col_nr=0;
        col_id=[];
        data=[];
        continue
			if col_nr<0 and len(tags)>0 and tags[0][0]=='_':
				#this is a field for the frame
				fkey=tags[0]
				if len(tags)>1:
					fval=tags[1]
					fields[fkey]=fval
				else: multi_line_field='START'
				continue
			if col_nr>=0 and len(tags)>0 and tags[0][0]=='_':
        col_id.append(tags[0]);
        col_nr+=1;
        continue
			if col_nr>=0 and len(tags)>0 and tags[0]=='stop_':
				loops.append( self.Loop( col_id, data ))
				col_nr=-1
        continue
			if col_nr>=0 and len(tags)>0:
        data.append(tags);
        continue
			if len(tags)>0 and tags[0]=='save_':
				return loops, fields
			if len(tags)>0 and tags[0][0]==';' and multi_line_field:
				if multi_line_field=='START':
					multi_line_field='CURRENT'
					mlf_data=[]
				elif multi_line_field=='CURRENT':
					multi_line_field=None
					fields[fkey]=mlf_data
				continue
			if len(tags)>0 and multi_line_field=='CURRENT':
				mlf_data.append(line[:-1])
		return loops, fields

	#go through all frames and store the respective loops
	def parse_file( self, file ):
		while True:
			SAVE_FRAME_CATEGORY, name=self.Next_Save_Frame( file );
			if SAVE_FRAME_CATEGORY=='NO_CATEGORY':
				break

			loops, fields = self.capture_loops( file )
			self.add_frame(name, loops, fields, SAVE_FRAME_CATEGORY )

	def add_frame(self, name, loops, fields, CATEGORY='GENERIC'):
		if not CATEGORY in self._frames:
			self._frames[CATEGORY]={}
		self._frames[CATEGORY][name]=self.Frame(name, fields, loops )


	#how many different frame categories have been found ?
	def nframes( self ):
		return len( self._frames )

	#extract columns according to _tag_dictionary from loop
	def _process_loop( self, loop, cols ):
		ind={}
		types={}
		#figure out which bmrb-type columns fit to the requested columns (in cols) in this loop
		for col in cols:
			for tag in _tag_dictionary[col].tags:
				if tag in loop.cols:
					ind[col]=loop.cols.index(tag);
					types[col]=_tag_dictionary[col].type
					break

		#should have now the indices of the requested columns and their type in ind and types
		#if no fitting columns return
		if len(ind)==0: return
#		print 'C', cols
#		print 'L',loop.cols
#		print 'I', ind
		#extract output dictionary
		output={}
		for col in cols:
			output[col]=[]

		#lines are already split into columns
		for line in loop.data:
			if line[0][0]=='#': continue
			for i,col in enumerate(cols):
#				print 'F', i, col, line
				if not col in ind.keys() or line[ind[col]]=='.':
					del cols[i]
					return self._process_loop( loop, cols )
				output[col].append(types[col](line[ind[col]]))
		return output

	#process frames of a certain category and return output that fits the requested columns
  def process_frame( self, category, cols ):
		for frame in self._frames[category].itervalues():
			outputs=[]
			for loop in frame.loops:
				output=self._process_loop( loop, copy.copy(cols) )
				if output:
					outputs.append((len(output),output))

		outputs=sorted(outputs,key=lambda x: x[0])
		return output


_SEQUENCE_FIELD_KEY='_Mol_residue_sequence'

def get_sequence( bmrb_file ):
	seq_frames=bmrb_file._frames['monomeric_polymer']
	sequences={}
	for frame in seq_frames.itervalues():
		try:
			sequences[frame.name]=''.join([ x.strip() for x in frame.fields[_SEQUENCE_FIELD_KEY]])
		except KeyError:
			raise  library.InconsistentInput("Cannot find field _Mol_residue_sequence_'. Maybe not a proper BMRB file?")
	return sequences[sequences.keys()[0]], sequences, len(sequences)

def read_cs_bmrb( filename, verbose=0 ):
	fasta=None
	bmrb_file = BmrbFile( filename )
	try:
		seq_frame=bmrb_file._frames['monomeric_polymer']
	except KeyError:
		seq_frame=None

	try:
		cs_frame=bmrb_file._frames['assigned_chemical_shifts']
	except KeyError:
		#reading under assumption of full STAR structure didn't work, reread from start this time assume that loops are at top-level
		loops, fields = bmrb_file.capture_loops( open(filename, 'r' ) )
	#	print 'L', loops
	#	print 'F', fields
		bmrb_file.add_frame('generic',loops,fields,'assigned_chemical_shifts');
		if not seq_frame and _SEQUENCE_FIELD_KEY in fields:
			fasta=''.join([ x.strip() for x in fields[_SEQUENCE_FIELD_KEY]])
			nmol=1
	if verbose and seq_frame:
		print "Molecule(s) found in BMRB:\n",seq_frame

	if not fasta:
		try:
			(fasta, fastas, nmol)=get_sequence( bmrb_file )
		except KeyError:
			fasta=None
			nmol=1

	if nmol>1:
		print "WARNING: found multiple molecules in BMRB"
		print fasta
		exit()

	return bmrb_file, fasta
