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
import argparse
import sys
import library
import noesy
import fasta
import traceback
#############################
#if len(argv) <=1:
#    Help()

#file = argv[1]
parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-rdc", help="a raw RDC table with A HN N RDC", required=True);
parser.add_argument("-tol", help="tolerances to match atom1 and atom2", nargs=2, type=float, default=[0.05, 0.5])
parser.add_argument("-atoms", help="atom types to match", nargs=2, default=['H','N'] )
parser.add_argument("-resonances", help="a file with resonances for assignment in prot format",required='True');
parser.add_argument("-fasta", help='sequence information, if not already contained in resonance file', default=None );
parser.add_argument("-target", help='provide target fasta sequence to work out trimming', default=None );
parser.add_argument("-o", help='output file with assigned peaks', default=None );
library.add_standard_args( parser )
args = parser.parse_args()

class RDC:
	class Ambiguous:
		pass
	class Unassigned:
		pass
	class Assignment:
		def __init__(self, rmatch, rself, tol ):
			self.r=rmatch
			self.rdc_freq=rself
			self.scores={}
			self.add_score( 'cs_match', -rmatch[0].pmatch(rself[0],tol[0])*rmatch[1].pmatch(rself[1],tol[1]) )

		def add_score(self, type, val ):
			self.scores[type]=val

		def total_score(self):
			return sum(self.scores.values())

		def resid(self):
			return self.r[0].resid()

		def freq(self,i):
			return self.r[i].freq()

		def delta(self,i):
			return self.r[i].freq()-self.rdc_freq[i]

		def __str__(self):
			return '(%3d, %5.3f %5.3f S=%5.3f)'%(self.resid(),self.delta(0), self.delta(1), self.total_score())

	def __init__(self,id,freq1,freq2,rdc):
		self.id=id
		self.freq=(freq1,freq2)
		self.rdc=rdc
		self._assignments=[]

	def __str__(self):
		s='RDC ( %s %5.3f %5.3f: %8s )  '%(self.id, self.freq[0],self.freq[1],self.rdc)
		s+=' '.join( [ '%s'%r for r in self._assignments ] )
		return s

	def append(self,r1,r2,tol):
		self._assignments.append(RDC.Assignment((r1,r2),self.freq,tol))

	def __len__(self):
		return len(self._assignments)

	def assignments(self):
		for i in self._assignments:
			yield i

	def restrict_to_single_assignment(self):
		if len( self._assignments ) > 1:
			sortass = sorted( self._assignments, key=RDC.Assignment.total_score )
			self._assignments=[sortass[0]]

	def raise_on_non_single(self):
		if len( self._assignments ) > 1:
			raise RDC.Ambiguous
		if len( self._assignments ) == 0:
			raise RDC.Unassigned

	def resid(self):
		self.raise_on_non_single()
		return self._assignments[0].resid()

	def total_score(self):
		self.raise_on_non_single()
		return self._assignments[0].total_score()

	def resonances(self):
		self.raise_on_non_single()
		return self._assignments[0].r

def dump_rdcs( rdc_by_resid ):
	for resid,rdclist in rdc_by_resid.iteritems():
		print resid, ' '.join(['%s'%rdc for rdc in rdclist ])

try:
#output:
	library.hello( __file__ )
	res_file = open( args.resonances, 'r')
	resonances=noesy.ResonanceList.read_from_stream( res_file )

	if args.fasta:
		sequence=fasta.read_fasta( args.fasta )

	if args.fasta and resonances.sequence():
		offset=fasta.find_fasta_offset( resonances.sequence(), sequence )
		if offset:
			raise library.InconsistentInput('Sequence in %s has offset of %d from sequence in %s'%(args.fasta, offset, args.resonances))
	if args.fasta:
		resonances.set_sequence(sequence)

	#read RDC raw data
	rdc_list=[]
	for line in open(args.rdc,'r'):
		tags=line.split()
		if len(tags)<1: continue
		if tags[0]=='ASS': continue
		if '#' == tags[0][0]: continue
		if len(tags)==4:
			try:
				rdc_list.append(RDC(tags[0],float(tags[1]),float(tags[2]),tags[3]))
#				print rdc_list[-1]
			except ValueError:
				print 'ignored line: ', line[:-1]

	#get initail matches --- anything within tolerances
	rdc_by_resid={}
	for rdc in rdc_list:
		for r1 in resonances.matches( args.atoms[0], rdc.freq[0], args.tol[0]):
			r2=resonances.by_atom( noesy.Atom( args.atoms[1], r1.resid() ) )
			if r2.match( rdc.freq[1], args.tol[1] ):
				rdc.append(r1,r2,args.tol)
				rdc_by_resid.setdefault(r1.resid(),[]).append(rdc)

	#resolve double assignments to a single RDC
	for resid,rdclist in rdc_by_resid.iteritems():
		rdc=rdclist[0]
		if len(rdclist)==1 and len(rdc)==1: continue
		if len(rdc)>1:
			for a in rdc.assignments():
				if a.resid()!=resid:
					a.add_score('other_places', 1.0/len(rdc_by_resid[a.resid()]))

	#dump
	dump_rdcs(rdc_by_resid)

	#rebuild rdc_by_resid and restrict to best-scoring assignment
	rdc_by_resid={}
	for rdc in rdc_list:
		rdc.restrict_to_single_assignment()
		try:
			rdc_by_resid.setdefault(rdc.resid(),[]).append(rdc)
		except RDC.Unassigned:
			print 'UNASSIGNED ',rdc
			pass

	#dump
	dump_rdcs(rdc_by_resid)

	#remove extra RDCs from residues that still have two RDCs assigned
	for resid,rdclist in rdc_by_resid.iteritems():
		if len(rdclist) > 1:
			#keep it a list but only one element
			s=sorted(rdclist, key=RDC.total_score )
			rdc_by_resid[resid]=[ s[0] ]
			if s[1].total_score()-s[0].total_score()<0.2:
				print 'difficult decision at residue %d: '%resid, ' '.join(['%s'%r for r in s])

	#dump
#	dump_rdcs(rdc_by_resid)

	if args.o:
		start=1
		end=99999
		if args.target:
			target_fasta=fasta.read_fasta(args.target)
			start=-fasta.find_fasta_offset( target_fasta,resonances.sequence() )+1
			end=start+len(target_fasta)-1
			sequence, end=fasta.cut_sequence(resonances.sequence(),start,end,1)

		file=open(args.o,'w')
		for resid,rdclist in rdc_by_resid.iteritems():
			rdc=rdclist[0]
			try:
				if resid>=start and resid<=end:
					file.write('%5d %3s  %5d %3s %8.3f\n'%(
							resid-start+1,
							rdc.resonances()[0].name(),
							resid-start+1,
							rdc.resonances()[1].name(),
							float(rdc.rdc)
							))
			except ValueError:
				print 'ignored assigned RDC %s'%rdc

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
