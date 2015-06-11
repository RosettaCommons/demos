#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-
from os.path import basename
import sys
import argparse

from os.path import exists
import traceback
#############################

parser = argparse.ArgumentParser(prog=basename(__file__), description="combines tempering stats from multiple replicas")
parser.add_argument("tempering", help="the tempering stats file written by rosetta canonical sampling module",default="tempering.stats");

args = parser.parse_args()

class Cell:
	def __init__(self, xchange, trials, accepts):
		self.xchange=xchange
		self.trials=trials
		self.accepts=accepts

	def combine( self, cell ):
		if self.xchange != cell.xchange:
			raise ValueError("Keyerature of cells incompatible, trying to combine key %s with key %s"%(self.xchange, cell.xchange))
		self.trials+=cell.trials
		self.accepts+=cell.accepts

	@classmethod
	def from_accept_ratio( obj, xchange, trials, accepts ):
		accepted_moves=int(round(1.0*accepts*trials))
		obj=Cell( xchange )
		obj.trials=trials
		obj.accepts=accepted_moves
		return obj

	def __str__( self ):
		if self.trials>0:
			return '%10s trials= %10d accepts= %8.4f'%(self.xchange,self.trials,1.0*self.accepts/self.trials)
		else:
			return '%10s trials= 0 accepts= nan'%(self.xchange)


file = open( args.tempering, 'r' )
levels={}
lines=file.readlines()
jobs={}
for line in lines:
	if line[0:5]=='Stats':
		tags=line.split()
		count=jobs.setdefault(tags[3],0)
		jobs[tags[3]]=count+1

counted_jobs=sorted( jobs.iteritems(), key=lambda x: x[1])

current_job=''
final_job=counted_jobs[-1][0]
for line in lines:
	tags=line.replace('=',' ').replace(';',' ').split()
	if tags[0]=='Stats':
		 current_job=tags[3]
	if current_job!=final_job: continue
	if tags[0]!='level': continue
	level=int(tags[1])
	temp=float(tags[3])
	xchange=tags[4]
	trials=int(tags[6])
	try:
		accept_ratio=float(tags[8])
		accepts=int(round(1.0*accept_ratio*trials))
	except ValueError:
		accepts=0

	#	cell=Cell.from_accept_ratio( xchange, trials, accepts )
	cell=Cell( xchange, trials, accepts )
	default_cell=Cell(xchange,0,0)
	default_cell.combine( cell )
	levels.setdefault(xchange,default_cell).combine( cell )

sort_keys=[]
for key,val in levels.iteritems():
	tags=key.split('_')
	skey=(int(tags[1]),int(tags[2]))
	if skey[0]<skey[1]:
		sort_keys.append( (skey,key,val) )

for sky,key,val in sorted( sort_keys, key=lambda x: x[0] ):
	back_key='HX_%d_%d'%(sky[1],sky[0])
	print '%s %s'%(val,levels[back_key])
