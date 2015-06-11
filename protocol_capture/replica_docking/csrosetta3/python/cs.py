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
import amino_acids
### toolbox library
#import library
from os import mkdir , makedirs

from warnings import *
import warnings

import traceback

from PDB import PDBParser
from PDB import PDBIO

from library import square

from library import mkdirp
import library

class NIH_table:
	def __init__(self):
		self.vars=None

		#format string
		self.format=""

		#line entries, indexed by (resi, atom)
		self.table={}

		#additional data entries:
		#for instance DATA SEQUENCE VAL would be data['SEQUENCE']='VAL'
		self.data={}

		self.has_atom=False

	def read_file( self, file):
		fd=open(file,'r')
		self.read( fd )
		return self

	def read( self, fd ):
		self.__init__()
		self.read_header( fd )
		if not self.vars or not self.format:
			raise MissingInput("cannot find VARS and FORMAT entries in file %s"%file )
		self.read_body( fd )
		return self

	def iteritems(self):
		return self.table.iteritems()

	def read_header( self, fd ):
		for line in fd:
			tags=line.split()
			if not tags: continue
			if "DATA"==tags[0]:
				key=tags[1];
				tags[1]=""
				if key in self.data.keys():
					tags[1]=self.data[key]
				self.data[key]="".join(tags[1:])
			if "VARS"==tags[0]:
				self.vars=tags[1:]
			if "FORMAT"==tags[0]:
			 	self.format=line[6:-1]
			if self.vars and self.format:
			 	break #finish reading-loop, do some more setup and resume reading file
		if not self.vars:
			raise library.MissingInput("cannot find VARS entry in header")
		if not self.format:
			raise library.MissingInput("cannot find FORMAT entry in header")
		if len(self.format.split())!=len(self.vars):
			raise library.InconsistentInput("VARS and FORMAT entry have different length: %s %s"%(self.vars, self.format))

	def read_body( self, fd ):
		floats=[]
		ints=[]
		f_tags=self.format.split()
		for f in f_tags:
			floats.append( "f" in f or "g" in f );
			ints.append( "d" in f );
		col_resid=self.vars.index("RESID")
		try:
			col_atom=self.vars.index("ATOMNAME")
			self.has_atom=True
		except:
			self.has_atom=False

		self.table={}
		for line in fd:
			tags=line.split()
			if len(tags)==len(self.vars):
				resi=int(tags[col_resid])
				try:
					atom=tags[col_atom]
					if atom=="HN":
						atom="H"
					key=(resi,atom)
				except:
					key=resi
				value=tags
				for ct in range(0, len(tags) ):
					if floats[ct]: value[ct]=float(tags[ct])
					elif ints[ct]: value[ct]=int(tags[ct])
				self.table[key]=value
			elif len(tags)>0:
				raise library.InconsistentInput("expected %d entries (from VARS) but found only %d in line %s.\n VARS: %s"%(len(self.vars),len(tags),line," ".join(self.vars)))

	def from_dict( self, dict ):
		self.vars=dict.keys()
		self.table={}
		for i,resi in enumerate(dict['RESID']):
			try:
				atom=dict['ATOMNAME'][i]
				self.has_atom=True
				key=(resi,atom)
			except:
				key=resi
			self.table[key]=[ tab[i] for tab in dict.itervalues() ]
		if 'RESNAME' in dict:
			fasta=get_sequence_from_column( self.get_slice('RESNAME') )
			self.data['SEQUENCE']=fasta
		return self

	def set_sequence( self, sequence ):
		sys.stderr.write('set sequence... \n')
		self.data['SEQUENCE']=sequence
		try:
			col_res1 = self.get_col('RESNAME')
		except ValueError:
			self.vars.append('RESNAME')
			self.format+=" %3s"
			col_res1 = None

		try:
			col_res3 = self.get_col('RESNAME3')
		except ValueError:
			self.vars.append('RESNAME3')
			self.format+=" %5s"
			col_res3 = None

		for key,val in self.table.iteritems():
#				val=self.table[key]
				try:
					resi=key[0]
				except:
					resi=key
				try:
					aa=sequence[resi-1]
				except IndexError:
					raise library.InconsistentInput('trying to find residue %d in fasta-string %s but string is too short'%(resi,sequence))
				aa3=amino_acids.short_to_long[aa]
				if col_res1: val[col_res1]=aa
				else: val.append(aa)
				if col_res3: val[col_res3]=aa3
				else: val.append(aa3)

		#renumber resid indices according to start...end
	def renumber( self, start, end ):
			col_resid=self.vars.index("RESID")
			new_table={}
			for key in self.table.keys():
				val=self.table[key]
				try:
					k=key[0]
					new_key=(key[0]-start+1, key[1])
					val[col_resid]=new_key[0]
				except:
					k=key
					new_key=key-start+1
					val[col_resid]=new_key
				if k<start or k>end:
					continue
				new_table[new_key]=val
			self.table=new_table

	def write_file( self, file ):
		fd=open(file,'w')
		self.write( fd )

	#write to file descriptor
	def write( self, fd ):
		self.write_header( fd )
		self.write_body( fd )

	def write_header( self, fd ):
		if not self.vars or not self.format: return
		for d in self.data.keys():
			out_str = self.data[d]
			for i in range(0,len(out_str),10):
				if i%50==0: fd.write("\nDATA %s"%d);
				fd.write(" "+out_str[i:(i+10)]);
			fd.write("\n");

		if self.data:
			fd.write("\n")
		fd.write("VARS %s\n"%(" ".join(self.vars)))
		fd.write("FORMAT %s\n\n"%(self.format))

	def write_body( self, fd ):
		keys=self.table.keys()
		keys.sort()

		for key in keys:
			entry=self.table[key]
			try:
				fd.write((self.format%tuple(entry))+"\n")
			except TypeError:
				print 'cannot write ', entry, ' with the given format ', self.format
				raise

	def get( self, resi, atom, key ):
		col=self.vars.index( key )
		return self.table[ (resi, atom) ][ col ]

	def put( self, resi, atom, key, value ):
		col=self.vars.index( key )
		self.table[ (resi, atom)][col]= value

	def get_col( self, key ):
		return self.vars.index( key )

	def has_col( self, key ):
		return key in self.vars

	def __contains__(self, key):
		return key in self.vars

	def copy( self, table_in ):
#		print self.vars
		col_map = [ table_in.get_col(v) for v in self.vars ]
		self.table={}
		for key in table_in.table.keys():
			self.table[key]=[ table_in.table[ key ][ col ] for col in col_map ]

	def get_slice( self, column ):
		col=self.get_col( column )
		slice={}
		for key in self.table.keys():
			slice[key]=self.table[key][col]
		return slice

	def add_slice( self, column, values, in_front_of=None, format='' ):
#		print column
#		print self.vars
		col=len(self.vars)
		insert=False
		if not column in self.vars:
			insert=True
			if in_front_of:
				col=self.get_col( in_front_of )
			self.vars.insert( col, column)
			if format:
				tags=self.format.split('%')
				tags.insert(col+1, format[1:]+' ' )
				self.format='%'.join( tags )
		else:
			col=self.get_col( column )
#		print col, self.vars

		for key in values.keys():
			if key in self.table:
				if insert:
					self.table[key].insert( col, values[key] )
				else:
					self.table[key][col]=values[key]
			else:
				fields=[0]*len(self.vars)
				colr=self.vars.index("RESID")
				colrn=None
				try:
					fields[colr]=key[0]
					cola=self.vars.index("ATOMNAME")
					fields[cola]=key[1]
				except:
					fields[colr]=key
				if "RESNAME" in self.vars:
					colrn=self.vars.index("RESNAME")
				if colrn and "SEQUENCE" in self.data:
					fields[colrn]=self.data["SEQUENCE"][fields[colr]-1]
				self.table[key]=fields
		return self

class TalosCSFile(NIH_table):
	def __init__(self):
		NIH_table.__init__(self)
		self.sequence=None
#		self=None

	def read_file(self, file ):
#		self=NIH_table()
		NIH_table.__init__(self)
		NIH_table.read_file( self, file )
		if 'SEQUENCE' in self.data:
			self.sequence=self.data['SEQUENCE']

	def write( self, fd ):
#lf:
#			return
		if self.sequence:
			self.data['SEQUENCE']=self.sequence
		NIH_table.write( self, fd )

	def write_file( self, file ):
		fd=open( file, 'w' )
		self.write(fd)

	def renumber( self, start, end ):
		self.sequence= self.sequence[(start-1):end]
		NIH_table.renumber( self, start, end )
		self.data['SEQUENCE']=self.sequence

	def from_table( self, table_in, sequence=None ):
#		self.table = NIH_table()
		NIH_table.__init__(self)
		if 'RESNAME' in table_in.vars:
			self.vars = ['RESID','RESNAME','ATOMNAME','SHIFT']
		else:
			raise library.MissingInput("require sequence information to make talos-table out of prot-table")
		self.format='%4d %1s %4s %8.3f'
		self.copy( table_in )
		self.sequence = sequence
		if not self.sequence and 'SEQUENCE' in table_in.data:
			self.sequence = table_in.data['SEQUENCE']
			self.data['SEQUENCE']=self.sequence
		#figure out if RESNAME needs translating from aa3-->aa1
		resn=self.get_slice( 'RESNAME' )
		if len(resn[resn.keys()[0]])==3:
			resn_translated=table_op( resn, resn, lambda x,y: amino_acids.longer_names[x.upper()] )
			self.add_slice( 'RESNAME', resn_translated )
		elif not len(resn[resn.keys()[0]])==1:
			raise library.InconsistentInput("RESNAME column in input table should have either aa3 or aa1 format, i.e., ALA or A for residue names")


def get_sequence_from_column( data ):
	min_resi = 1000000
	max_resi = -min_resi
	for key,a in data.iteritems():
		if min_resi > key[0]: min_resi = key[0]
		if max_resi < key[0]: max_resi = key[0]
	sequence=list("-"*(max_resi))
	for key,a in data.iteritems():
		if len(a)==1:
			sequence[key[0]-1]=a
		else:
			sequence[key[0]-1]=amino_acids.short2long( a )
	return ''.join(sequence)

class ProtCSFile(NIH_table):
	def __init__(self):
		NIH_table.__init__(self)
		self.sequence=None
#		self.table=None

	def read(self, fd, has_aa3, has_aa1, sequence=None, header=False ):
		#self.table=NIH_table()
		NIH_table.__init__(self)
		if not header:
			self.vars=['INDEX','SHIFT','SIGMA','ATOMNAME','RESID']
			self.format='%8d %10.3f %10.3f %7s %8d'

			if has_aa3:
				self.vars.append('RESNAME3')
				self.format=self.format+' %6s'
			if has_aa1:
				self.vars.append('RESNAME')
				self.format=self.format+' %4s'
		else: #has header
			self.read_header( fd )
			has_aa3 = 'RESNAME3' in self.vars
			has_aa1 = 'RESNAME' in self.vars

		#read body
		self.read_body( fd )

		if not has_aa1:
			if has_aa3:
				aa3=self.get_slice('RESNAME3')
				aa={}
				for key,a in aa3.iteritems():
					aa[key]=amino_acids.longer_names[ a.upper() ]
				self.add_slice('RESNAME', aa)
				self.format=self.format+" %5s"
			elif sequence:
				aa3=self.get_slice('INDEX')
				aa={}
				for key,a in aa3.iteritems():
					aa[key]=sequence[key[0]-1]
				self.add_slice('RESNAME', aa)
				self.format=self.format+" %5s"
		if not sequence and (has_aa1 or has_aa3):
			aa1=self.get_slice('RESNAME')
			sequence = get_sequence_from_column( aa1 )

		if 'SEQUENCE' in self.data:
			sequence = self.data['SEQUENCE']
		return sequence

	def renumber( self, start, end ):
		if self.sequence:
			self.sequence= self.sequence[(start-1):end]
			self.data['SEQUENCE']=self.sequence
		NIH_table.renumber( self, start, end )



  def read_file( self, file, sequence=None ):
		fd=open(file,'r')
		self.read_stream( fd, sequence )

	def read_stream( self, fd, sequence=None ):
		has_header = False
		has_aa1 = False
		has_aa3 = False
		for l in fd:
			tags=l.split()
			if len(tags)<1: continue
			if 'VARS' in tags[0]:
				has_header = True
				break

#			print len(tags), has_aa3, has_aa1, len(tags[5]), len(tags[len(tags)-1])
			if len(tags)==6 and len(tags[5])==3:
				has_aa3 = True
			if len(tags)>=6 and len(tags[len(tags)-1])==1:
				has_aa1 = True
			if len(tags)>6 and len(tags[len(tags)-2])==3:
				has_aa3 = True

	#	assert( fd.seekable() )
		fd.seek(0)
#		print has_header, has_aa1, has_aa3
		self.sequence = self.read(fd, has_aa3, has_aa1, sequence=sequence, header=has_header )



	def from_table( self, table_in, sequence=None ):
		NIH_table.__init__(self)
		if 'RESNAME' in table_in.vars:
			self.vars = ['SHIFT','ATOMNAME','RESID','RESNAME']
		self.format='%10.3f %7s %8d %4s'
		self.copy( table_in )
		self.sequence = sequence

		lines=self.get_slice( 'RESID' )
		indices={}
		sigmas={}
		keys=lines.keys()
		keys.sort()
		for id,key in enumerate(keys):
			indices[key]=id+1
			sigmas[key]=0.00
		self.add_slice( 'INDEX', indices, in_front_of='SHIFT', format='%8d' )
		self.add_slice( 'SIGMA', sigmas, in_front_of='ATOMNAME', format='%5.3f' )

		if not self.sequence and 'SEQUENCE' in table_in.data:
			self.sequence = table_in.data['SEQUENCE']
			self.data['SEQUENCE']=self.sequence


	def write( self, fd, header=False ):
		if header:
			if self.sequence:
				self.data['SEQUENCE']=self.sequence
			self.write_header( fd )
		self.write_body( fd )


#SpartaCalculator: allows to spawn Sparta calculations:
#	  refCS as reference chemical shifts
#   pdb the input pdb structure
#   pred_dir a directory name to store all Sparta data files, these are automatically reused
#   for subsequent runs with same pdb structure, unless overwrite=True.
#   if the Sparta class is instantiated for same pdb multiple times it will be copied from
#   previous instance using the cached 'SpartaStore'

class SpartaCalculator:
	def __init__( self, refCS, pred_dir, overwrite=False ):
		self.refCS=refCS
		self.pred_dir=pred_dir
		self.overwrite=overwrite
		self.SpartaStore={}
	def get_sparta( self, pdb ):
		return self.Sparta( pdb, self.refCS, self.pred_dir, self.overwrite, self.SpartaStore )

	class Sparta:
	  def __init__( self, pdb, refCS, pred_dir, overwrite, SpartaStore ):
		  if pdb in SpartaStore:
				old=SpartaStore[pdb]
#				print "retrieve: ", pdb, id(old)
				self.__copy_constr__( old )
				self.refCS = refCS
			else:
				library.mkdirp( pred_dir )
				self.refCS = refCS
				self.pred_dir = pred_dir
				self.pdb = pdb
				self.name= splitext( basename(pdb) )[0]
				self.overwrite = overwrite
				self.pred_fn = join( self.pred_dir, "pred_"+self.name+".tab" )
				self.struct = join( self.pred_dir, "struct_"+self.name+".tab" )
				self.csCalc_fn = join( self.pred_dir, "csCalc_"+self.name+".tab" )

				self.pred = NIH_table()
				self.csCalc = NIH_table()

				if not self.results_exist() or self.overwrite:
					self.runSparta()
				self.read_results()
				SpartaStore[pdb]=self;
			#print "generate: ", pdb, id(self)


		def __copy_constr__( self, old ):
			self.refCS = old.refCS
			self.pred_dir = old.pred_dir
			self.pdb = old.pdb
			self.name = old.name
			self.overwrite = old.overwrite
			self.pred_fn = old.pred_fn
			self.struct = old.struct
			self.csCalc_fn = old.csCalc_fn
			self.pred = old.pred
			self.csCalc = old.csCalc

		def __str__( self ):
			return "SPARTA: refCS: %s, pred_dir: %s, pdb: %s\n pred: %s, struct: %s, csCalc: %s\n"%(self.refCS,\
                                                         self.pred_dir, self.pdb, self.pred_fn, self.struct, self.csCalc_fn )

		def runSparta( self ):
			subprocess.call(["sparta+","-in %s -out %s -outS %s -ref %s -offset -outCS %s"%( self.pdb, self.pred_fn, self.struct, self.refCS, self.csCalc_fn ) ], stderr=subprocess.PIPE);

		def read_results( self ):
			if exists( self.pred_fn ):
				self.pred.read_file( self.pred_fn )
			if exists( self.csCalc_fn ):
				self.csCalc.read_file( self.csCalc_fn )

		def results_exist( self ):
			return exists( self.pred_fn ) and ( exists( self.csCalc_fn ) or not self.refCS )

		def get_slice( self, column ):
			try:
				return self.pred.get_slice( column )
			except ValueError:
				try:
					return self.csCalc.get_slice( column )
				except ValueError:
					raise library.MissingInput("field %s is neither found in pred_ nor in csCalc_ files of Sparta+ output"%(column) )

#End Sparta class definition


def table_op( table1, table2, func, selection=None, default=0 ):
    table={}
    for key in table1.keys():
        try:
					if not selection or ( key[0] in selection.resi ):
            table[key]=func( table1[key], table2[key] )
					else:
						table[key]=default
        except KeyError:
            pass
    return table


def average_shifts( column ):
    global args
    pdblist=args.pdbs;
    pdb = pdblist[ 0 ]
    sparta = Sparta( pdb )
    tab=sparta.get_slice( column );
    #setup empty tables
    av_shifts=table_op( tab, tab, lambda x,y: 0.0 );
    av2_shifts=av_shifts;
    inv_N = 1.0 / len(pdblist)
    for pdb in pdblist:
        sparta = Sparta( pdb )
#        print "got : ", id( sparta )
        new_tab=sparta.get_slice( column )
        av_shifts=table_op( av_shifts, new_tab, lambda x,y: x+y*inv_N )
        av2_shifts=table_op( av2_shifts, new_tab, lambda x,y: x+y*y*inv_N )
    #subtract square of average : std = <x^2>-<x>^2
    std_shifts=table_op( av2_shifts, av_shifts, lambda x,y: x-y*y )
    return av_shifts,std_shifts

#by default no extra weighting per atom, since Sparta+ results, for instance, are already weighted according to sigma
def write_to_bfactor( pdbin, pdbout, result, sum_residue=True, \
                      func=(lambda x: x), atoms={"N":1.0,"H":1.0,"C":1.0,"CB":1.0,"CA":1.0,"HA":1.0}, \
                      response=(lambda x: x)
                      ):

	with warnings.catch_warnings(record=True) as w:
		parser=PDBParser()
		s=parser.get_structure( "t000_", pdbin )

	final_values={}
	for chain in s[0]:
		for residue in chain:
			resid=residue.get_id()[1]
			sum=0;
			if sum_residue:
				for atom in atoms:
					try:
						sum=sum+func(result[(resid,atom)])*atoms[atom]
					except KeyError:
						pass
				for atom in residue:
					atom.set_bfactor( response( sum ) )
				final_values[resid]=response( sum )
			else:
				for atom in residue:
					try:
						v=response( result[(resid,atom.get_name())] )
						atom.set_bfactor( v )
						final_values[atom]=v
					except KeyError:
						atom.set_bfactor( 0.0 )



	io=PDBIO()
	io.set_structure( s )
	io.save( pdbout )
	return final_values

def sum_table( table, selection ):
	sum=0
	for key in table.keys():
		if key[0] in selection.resi:
			sum = sum + table[key]
	return sum

class AtomSelection:
	def __init__( self, selection ):
		self.selection=selection
		self.resi=[]
		if selection:
			tagsp=selection.split('+')
			for t in tagsp:
				if '-' in t:
					tags=t.split('-')
					for i in range(int(tags[0]),int(tags[1])+1):
						self.resi.append(i)
				else:
					self.resi.append(int(t))
		else:
			self.resi=range(1, 1000)

