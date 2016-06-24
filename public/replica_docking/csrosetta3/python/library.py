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

from os.path import basename
from os import mkdir , makedirs
import string
import amino_acids
import textwrap
import hashlib


def square(x):
	return x*x

def hashfile(afile, hasher, blocksize=65536):
	buf = afile.read(blocksize)
	while len(buf) > 0:
		hasher.update(buf)
		buf = afile.read(blocksize)
	return hasher.digest()

def diff_files_md5(file1,file2):
	return hashfile(open(file1,'r'),hashlib.md5() ) == hashfile(open(file2,'r'),hashlib.md5())


def hello( name ):
	print '-'*75, """
---------             CS-Rosetta 3.0   (toolbox)                -----------
-----                                                                 -----
--     website: www.csrosetta.org                                        --
--     copyright: Oliver Lange                                           --
--     Reference: Lange et al. PNAS 2012,                                --
--         www.pnas.org/cgi/doi/10.1073/pnas.1203013109                  --
--                                                                       --
--     program name: %30s                      --
--                                                                       --
-----                                                                 -----
--------                                                            -------"""%basename(name)
	print '-'*75

def add_standard_args( parser ):
	parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )

class Tracer:
	def __init__(self, channel):
		self.channel = channel+": "

	def out( self, msg ):
		pass
#		print self.channel+msg


class LibException(Exception):

	def __init__(self,msg):
		Exception.__init__(self)
		self.msg=[msg]

	def __str__(self):
		#        return "\n------------------      %s      --------------------------------------------\n"%(self.msg)
		str = ""
		for msg in self.msg:
			str = str + "\n>> "+"\n>> ".join( textwrap.wrap( msg ) )
		return str

	def add( self, str ):
		self.msg.append( str )

class StubbedOut(LibException):
	pass

class MissingInput(LibException):
	pass

class InconsistentInput(LibException):
	pass

class RunException(LibException):
	pass

class MethodException(LibException):
	def __init__(self,method,msg):
		LibException.__init__(self,msg)
		self.method=method
	def __str__(self):
		return "\n\n"+"-"*80+"\n"\
			 +"-------------------           in Method: "+self.method.__str__()+"          -----------------------\n"\
			 +self.method.module_description()+"\n"\
			 +"-"*80+"\n"\
			 +LibException.__str__(self)+"\n\n"\
			 +"-"*80+"\n"\

class ProgramError(LibException):
	def __str__(self):
		return """
\n\n------------------      %s      --------------------------------------------\n
\n please send a description of what you were doing, the cmd-line and the stack-trace (run with -traceback) to
 to : oliver.lange@tum.de \n
"""%(self.msg)

#returns the full name of prog prog.default.linuxgccrelease
# takes the fist present, limit using the option extras if you only like mpi for instance
def rosetta_executable( prog, extras=['default','static','mpi'] ):
	import subprocess
	platforms=['linux','macox','windows']
	cc=['gcc','icc']
	debug=['release','debug']
	testcmd='type %s >/dev/null 2>&1'
	for di in debug:
		for ci in cc:
			for ei in extras:
				for pi in platforms:
					name='%(prog)s.%(ei)s.%(pi)s%(ci)s%(di)s'%locals()
					try:
						subprocess.check_call(testcmd%name,shell=True)
						return name
					except subprocess.CalledProcessError as exc:
						pass

def read_rigid_file(rigid_file):
	lines=open(rigid_file,'r').readlines()
	start=100000
	end=0
	for line in lines:
		tags=line.split()
		if not len(tags): continue
		if tags[0]!='RIGID': raise InconsistentInput('expected RIGID at start of line %s in file %s'%(line,rigid_file))
		s=int(tags[1])
		e=int(tags[2])
		if start>s: start=s
		if end<e: end=e
	return start,end

def obj_is_list( val ):
	is_list=not isinstance(val,basestring)
	try:
		for file in val:
			pass
	except:
		is_list=False
	return is_list



def cst_is_centroid(file):
	lines = open(file).readlines()
	for line in lines:
		if (line[0]=="#"): continue;
		l=string.split(line);
		try:
			resi1 = int( l[2] );
			resi2 = int( l[4] );
		except:
			print '\n\n expected integer in col 3 and 5 of:\n'
			print l[2],l[4],'\n',l
			exit()

		atom1=l[1];
		atom2=l[3];
		centroid=["H","HA","CA","CB","C","O","N","ZN"]
		if not atom1 in centroid: return False
		if not atom2 in centroid: return False
	return True

def read_aa3_sequence(file):
	lines=open(file,'r').readlines()
	sequence=""
	for line in lines:
		sequence=sequence+amino_acids.longer_names[ line[:-1] ]
	return sequence

def mkdirp( target_dir ):
	try:
#		if ".." in target_dir:
#			mkdir( target_dir );
#		else:
		makedirs( target_dir );
	except OSError as err:
		if not err[1]=="File exists":
			print err;
			raise

def _flatten_list(val, list ):
	if isinstance(val,basestring):
		list.append(val)
		return
	try:
		for i in val:
			_flatten_list(i, list)
	except:
		list.append(val)

def flatten_list( val ):
	list=[]
	_flatten_list( val, list)
	return list

def is_gzipped( file ):
	return file[-3:]==".gz"

def cst_has_HA_atoms( file ):
	lines = open(file).readlines()
	for line in lines:
      if (line[0]=="#"): continue
      l=string.split(line);
      try:
			resi1 = int( l[2] );
			resi2 = int( l[4] );
      except:
			print '\n\n expected integer in col 3 and 5 of:\n'
			print l[2],l[4],'\n',l
			exit()

		atom1=l[1];
      atom2=l[3];
      ha_atom=["HA"]
      if atom1 in ha_atom: return True
      if atom2 in ha_atom: return True
	return False

def _upl2mini_line( line, sequence_offset, max_res, sep=4, pad=0.15 ):
	if line[0]=="#": return ""
	cols=line.split()
	if (len(cols)<5): return ""
	resi=int(cols[0])
	resj=int(cols[3])
	atomi=cols[2]
	atomj=cols[5]
	atomi=atomi.replace("HN","H").replace("M","Q")
	atomj=atomj.replace("HN","H").replace("M","Q")
	if  ( (resi-resj)>=sep or (resj-resi)>=sep ) and (resi-sequence_offset)>0 and (resj-sequence_offset)>0 and (resi-sequence_offset<=max_res) and (resj-sequence_offset <= max_res):
		return "AmbiguousNMRDistance %5s %5d %5s %5d BOUNDED 1.5 %5.3f 0.3 NOE; rawdata %5.3f\n"%(atomi,int(cols[0])-sequence_offset,atomj,int(cols[3])-sequence_offset,float(cols[6])+pad,float(cols[6]))
	else: return ""


# def upl2mini_find_offset( fasta, lines ):
# 	#work out offset first

# 	#build dictionary res-number --> aa3
# 	dict={}
# 	for line in lines:
# 		if line[0]=="#": continue
# 		cols=line.split()
# 		if len(cols)<5: continue
# 		try:
# 			cols[1]
# 			cols[1][0:3]
# 		except:
# 			print 'Exception in line: %s'%line, cols

# 		if int(cols[0]) not in dict:
# 			dict[int(cols[0])]=cols[1][0:3]
# 		else:
# 			if dict[int(cols[0])]!=cols[1][0:3]:
# 				raise InconsistentInput("residue names within the upl file are inconsistent -- maybe to files with different sequence offset have been merged ? ")

# 		if int(cols[3]) not in dict:
# 			dict[int(cols[3])]=cols[4][0:3]
# 		else:
# 			if dict[int(cols[3])]!=cols[4][0:3]:
# 				raise InconsistentInput("residue names within the upl file are inconsistent -- maybe to files with different sequence offset have been merged ? ")



# 	#residue numbers as "keys"
# 	keys=dict.keys()
# 	keys.sort()
# 	print "highest residue number", keys[-1] #last residue
# 	max_upl_res=keys[-1]

# 	#create fasta sequence: use - for positions where aa3 is not known
# 	upl_fasta=list("-"*keys[-1])
# 	for key in keys:
# 		try:
# 			upl_fasta[key-1]=amino_acids.longer_names[ dict[key] ]
# 		except:
# 			print "WARNING: could not understand residue name %s"%dict[key]
# 	upl_fasta="".join(upl_fasta)
# 	print "TARGET FASTA: %s\n   UPL FASTA: %s\n"%(fasta, upl_fasta)

# 	#find offset
# 	start_upl=keys[0]-1
# 	inconsistent=False

# 	#find first usable pieces of sequence
# 	last_gap=0
# 	first_aa=0
# 	matchable_pieces={}
# 	for res_upl in range( 0, len(upl_fasta) ):
# 		#if this is a gap --> then is the end of the next matchable_piece
# 		if upl_fasta[res_upl]=='-':
# 			if last_gap<first_aa and res_upl-first_aa>4:
# 				matchable_pieces[upl_fasta[first_aa:res_upl]]=first_aa
# 			last_gap=res_upl
# 		elif last_gap == res_upl-1:
# 			first_aa=res_upl

# 	#start with longer pieces first to find match for offset
# 	sorted_pieces=sorted(matchable_pieces,key=len,reverse=True)

# 	#go through matchable pieces and see if it can be found in fasta
# 	inconsistent=True
# 	sequence_offset=0
# 	for pi in sorted_pieces:
# 		if pi in fasta:
# 			sequence_offset=matchable_pieces[pi]-fasta.index( pi )
# 			inconsistent=False
# 			break

# 	if len(sorted_pieces):
# 		print pi, " is in ", fasta
# 		print "with offset ", sequence_offset
# 		#now go from back to fron along the original fasta sequence to look for a matching piece
# 		if not inconsistent:
# 			print "found offset: %d\n"%sequence_offset
# 		else:
# 			print "cannot determine offset\n"
# 			raise InconsistentInput("upl fasta is inconsistent with the target fasta sequence")
# 	else:
# 		print "not sufficient sequence information in upl file -- assuming no offset"

# 	start_upl=sequence_offset-1
# 	#final-check:
# 	#all residue positions known from upl-fasta against the given fasta sequence using our offset
# 	inconsistent=False
# 	for ct in range( 0, min(len(fasta), max_upl_res-sequence_offset) ):
# 		start_upl=start_upl+1
# 		if upl_fasta[start_upl]=='-':
# 			continue
# 		if not fasta[ct]==upl_fasta[start_upl]:
# 			inconsistent=True
# 			print "found inconsistent aa %c-%c at position %d"%(fasta[ct],upl_fasta[start_upl],ct+1)


# 	if inconsistent:
# 		raise InconsistentInput("upl fasta is inconsistent with the target fasta sequence")

# 	return sequence_offset

def upl2mini( file, output_fn, fasta=None, sequence_offset=0, QFall_file=None, sep=4, pad=0.15, verbose=1 ):
	lines=open(file).readlines()
	length=200000
	from fasta import upl2fasta, find_fasta_offset
	if fasta:
#		sequence_offset=upl2mini_find_offset( fasta, lines )
		upl_fasta=upl2fasta( lines )
		sequence_offset=find_fasta_offset( upl_fasta, fasta, verbose=verbose )
		length=len(fasta)

	QFall=[]
	noQF=[]
	for line in lines:
		if "#QF" in line:
			QFall.append(line)
		else:
			noQF.append(line)

	if len(noQF) > 0 or not QFall_file:
		if isinstance(output_fn, str):
			output=open(output_fn,'w')
		else:
			output=output_fn
		for line in noQF:
			output.write(_upl2mini_line( line, sequence_offset, length, sep, pad  ))


   if len(QFall) > 0 and QFall_file:
		output=open(QFall_file,'w')
		for line in QFall:
			output.write(_upl2mini_line( line, sequence_offset, length, sep, pad ))

	return len(noQF),len(QFall)

