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

from os.path import basename
import os
import string
import argparse
import sys
#sys.path.append('/home/zak/Downloads/numpy-1.6.2')
import library
import StringIO
import noesy
#from noesy import CrossPeakInfo
#from noesy import CrossPeakList
#from noesy import Resonance
#from noesy import ResonanceList
from noesy import SpinSystem, get_strips,spinsystem_evaluation
import fasta


parser =argparse.ArgumentParser(description="assign noesy peaks",
                                add_help=True)
parser.add_argument("-prot", help="chemical shift file",default=None);
parser.add_argument("-peaks",nargs='*', help='files with peak-lists in xeasy format',default=None);
parser.add_argument("-resid",help='whichr residue is validated',type=int, default=1);
#parser.add_argument("-putative_num", help='how many putative stripes will be generated for one proton',type=int, default=50 )
parser.add_argument("-new_label", help='Besides HA-CA and HB-CB, another label atom', default='CG' )
parser.add_argument("-new_H", help='Besides HA-CA and HB-CB, another label atom', default='HG2' )
parser.add_argument("-lib", help="library of CS distribution bounds",default=None )
parser.add_argument("-fasta", help="fasta file",default=None )
parser.add_argument("-score", help="score",type=float,default=0.5 )
parser.add_argument("-score_type", help="score type",default='simple' )
library.add_standard_args( parser )
args = parser.parse_args()

#librarylist=open(args.library,'r').readlines()
if 'csrosettaDir' in os.environ:
	library_file=os.environ['csrosettaDir']+"/database/cs_distribution.txt"
else:
	library_file=args.lib
resonance_list=noesy.ResonanceList.read_from_stream( open(args.prot,'r') )
crosspeaks=noesy.read_peak_files(args.peaks)

lib=open(library_file,'r').readlines()
if args.fasta:
	seq=fasta.read_fasta(args.fasta)

systems=[]
strips=get_strips.get_strips(args.resid,seq,lib,resonance_list,crosspeaks)
# for HA in HAs:
# 	for HB in HBs:
# 		sp=SpinSystem(args.resid)
# 		sp.add_strip(HN)
# 		sp.add_strip(HA)
# 		sp.add_strip(HB)
# 		if sp.score(args.score_type) > args.score:
# 			systems.append(sp)

for r in strips:
	sp=SpinSystem(args.resid)
	#print "length of r",len(r)
	for g in r:
		sp.add_strip(g)
	if sp.score(args.score,args.new_H)>args.score:
		systems.append(sp)
		prob=spinsystem_evaluation.spinsystem_evaluation(sp,lib,seq[args.resid-1])
		print "the incorrect prob of spinsystem is %(1)5.3e and the original score is %(2)5.3e\n"%{'1':prob,'2':sp.score(args.score,args.new_H)}


# systems=sorted(systems, key=SpinSystem, reverse=True )
# for sys in systems:
# 	print '%5.3f %s'%(sys.score(args.score,args.new_H), sys)
