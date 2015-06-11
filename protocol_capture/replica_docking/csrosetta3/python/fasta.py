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
import string
import amino_acids
import library
import sys
def cut_sequence(sequence,start,end,verbose):
	if ( end==0 ):
		end=len(sequence)
	if ( end > len(sequence) ):
		sys.stderr.write("WARNING: your end value is out of range\n");
	new_seq=sequence[(start-1):end];
	if verbose:
		print "\nold fasta:\n %s"%sequence
		print "\nnew fasta:\n %s%s%s"%('-'*(start-1),new_seq,'-'*(len(sequence)-end))
	return new_seq, end

def read_fasta(file):
	lines=open(file,'r').readlines()
	for line in lines:
		if line[0]==">":
			continue
	return line[:-1];

def write_fasta(file, fasta, tag=None):
	fd=open(file,'w')
	if not tag:
		tag='t000_'
	fd.write(">%s\n"%tag)
	fd.write("%s\n"%fasta)
	fd.close()

def find_fasta_offset( fasta1, fasta2, verbose=1 ):
	try:
		return _find_fasta_offset( fasta1, fasta2, verbose )
	except library.InconsistentInput:
		return -_find_fasta_offset( fasta2, fasta1, verbose )

def _find_fasta_offset( fasta1, fasta2, verbose=1 ):
	if verbose>0:
		print "trying to match the following sequences"
		print fasta1
		print fasta2

	#find first usable pieces of sequenc
	last_gap=0
	first_aa=0
	matchable_pieces={}
	fasta1=fasta1+'-'
	for res_upl in range( 0, len(fasta1) ):
		#if this is a gap --> then is the end of the next matchable_piece
		if fasta1[res_upl]=='-':
			if last_gap<first_aa and res_upl-first_aa>4:
				matchable_pieces[fasta1[first_aa:res_upl]]=first_aa
			last_gap=res_upl
		elif last_gap == res_upl-1:
			first_aa=res_upl

	#start with longer pieces first to find match for offset
	sorted_pieces=sorted(matchable_pieces,key=len,reverse=True)

	#go through matchable pieces and see if it can be found in fasta
	inconsistent=True
	sequence_offset=0
	for pi in sorted_pieces:
		for l in range( 5, len(pi)):
#			print ' try to find '+pi[0:l]+' in ',fasta2
			if pi[0:l] in fasta2:
				sequence_offset=matchable_pieces[pi]-fasta2.index( pi[0:l] )
				inconsistent=False
				break

	if len(sorted_pieces):
		if verbose>0: print pi, " is in ", fasta2
		if verbose>0: print "with offset ", sequence_offset
		#now go from back to fron along the original fasta2 sequence to look for a matching piece
		#if not inconsistent:
		#	print "found offset: %d\n"%sequence_offset
		if inconsistent:
			print "cannot determine offset\n"
			raise library.InconsistentInput("source fasta is inconsistent with the target fasta sequence")
	else:
		raise library.MissingInput("not sufficient information to determine offset")
	return sequence_offset

def upl2fasta( lines ):
  #work out offset first
	dict={}
	for line in lines:
		if line[0]=="#": continue
		cols=line.split()
		if len(cols)==0: continue
		dict[int(cols[0])]=cols[1][0:3]
		dict[int(cols[3])]=cols[4][0:3]

	#residue numbers as "keys"
	keys=dict.keys()
	keys.sort()
#print "highest residue number", keys[-1] #last residue
	max_upl_res=keys[-1]

	#create fasta sequence where AA is known -- otherwise
	upl_fasta=list("-"*keys[-1])
	for key in keys:
		upl_fasta[key-1]=amino_acids.longer_names[ dict[key] ]

	upl_fasta="".join(upl_fasta)
	return upl_fasta

def pdb2fasta( file ):
	lines=open(file,'r').readlines()
	fasta=[];
	for line in lines:
		ss=string.split( line );
		if ss[0]!='ATOM': continue
		if len(ss)>5 and ss[2]=="CA":
			fasta.append(amino_acids.longer_names[ss[3]]);
	return "".join(fasta)
