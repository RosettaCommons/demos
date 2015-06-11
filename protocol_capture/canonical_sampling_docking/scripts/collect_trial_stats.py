#!/usr/bin/env python2.7
## -*- mode:python;tab-width:2;indent-tabs-mode:t;rm-trailing-spaces:t;python-indent:2 -*-
from os.path import basename
import sys
import argparse

from os.path import exists
import traceback
#############################

parser = argparse.ArgumentParser(prog=basename(__file__), description="combines tempering stats from multiple replicas")
parser.add_argument("tempering", help="the tempering stats file written by rosetta canonical sampling module",default="trial.stats");

args = parser.parse_args()

file = open( args.tempering, 'r' )

levels={}
ct_levels={}
class Cell:
	def __init__(self, move, temp, trials, accepts):
		self.move=move
		self.temp=temp
		self.trials=trials
		self.accepts=accepts

	def combine( self, cell ):
		if self.move != cell.move or self.temp != cell.temp:
			raise ValueError("Keyerature of cells incompatible, trying to combine key %s with key %s"%(self.move, cell.move))
		self.trials+=cell.trials
		self.accepts+=cell.accepts

	def __str__( self ):
		if self.trials>0:
			return '%20s temp= %8.3f  trials= %10d accepts= %8.4f'%(self.move, self.temp, self.trials,1.0*self.accepts/self.trials)
		else:
			return '%10s temp= %8.3f  trials= 0 accepts= nan'%(self.move, self.temp)

for line in file:
	tags=line.replace('=',' ').replace(';',' ').split()
	if tags[0]!='level': continue
	level=int(tags[1])
	temp=float(tags[3])
	move=tags[4]
	trials=int(tags[6])
	try:
		accept_ratio=float(tags[8])
		accepts=int(round(1.0*accept_ratio*trials))
	except ValueError:
		accepts=0

	#	cell=Cell.from_accept_ratio( move, trials, accepts )
	cell=Cell( move, temp, trials, accepts )
	default_cell=Cell(move,temp,0,0)
	default_cell.combine( cell )
	levels.setdefault(move+"_%d"%level,default_cell).combine( cell )
	sum_trial = ct_levels.setdefault(level, 0 )
	sum_trial+=trials
	ct_levels[ level ] = sum_trial

for key in sorted( levels.keys() ):
	print levels[key]
for key, value in ct_levels.iteritems():
	print 'level= %5d trials= %10d'%(key, value)
