#!/usr/bin/env python2.7
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

import string
from silent_lib import ReadSilentData
import sys

from os.path import basename
import argparse
from ExampleArgumentParser import ExampleArgumentParser
import traceback
import library

from sys import argv,exit
from os import popen, system
from os.path import basename
import string

parser = ExampleArgumentParser(prog=basename(__file__), description="file silent file to extract decoys with lowest/highest score values",
examples=['%(prog)s decoys.out -score 10 > low_10.out',
          '%(prog)s decoys.out -rms 10 > low_10_by_rms.out'],
aliases=['extract_decoys'])

parser.add_argument("silentfile", nargs='*', help="input silent file", default=None);
parser.add_argument('columns', nargs='*', help='columns for output');
parser.add_argument('-formula', help='give a formula including score-column names: score+max(300,chem_shifts)', default=None )
parser.add_argument('-R',help='reverse output', action='store_true', default=False)
parser.add_argument('-N',help='number of decoys to write', type=int, default=None)

args = parser.parse_known_args()[0]
from sys import argv,exit
from os import popen, system
from os.path import basename
import string
import subprocess

REVERSE=''
NSTRUCT=10
if not args.formula:
	try:
    NSTRUCT = int(argv[-1])
    del(argv[-1])
	except:
    NSTRUCT = 2

	scorecol_defined = 0
	try:
    SCORECOL = int(argv[-1])
    del(argv[-1])
    scorecol_defined = 1
	except:
    SCORECOL = -1

	REVERSE = ''
	if SCORECOL > 0:
    REVERSE = ' --reverse '

#Another possibility... user supplies -rms or +rms
	scorecol_name_defined = 0
	if not scorecol_defined:
    scorecol_name = argv[-1]
    if scorecol_name[0] == '-':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        del( argv[-1] )
        REVERSE = ''
    if scorecol_name[0] == '+':
        scorecol_name_defined = 1
        scorecol_name = scorecol_name[1:]
        REVERSE = '-r'
        del( argv[-1] )

infiles = args.silentfile
if args.R:
	if REVERSE=='':
		REVERSE='-r'
	else:
		REVERSE=''

if args.N:
	NSTRUCT=args.N

#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

for infile in infiles:
    tags = []
    scoretags = string.split( popen('head -n 2 '+infile).readlines()[1] )
		if args.formula:
			used_tags=[ tag for tag in scoretags if tag in args.formula ]
			sys.stderr.write("used columns: %s\n"%" ".join(used_tags))
			if len(used_tags)==0:
				raise library.InconsistentInput("no valid column name found in %s\nChoose from %s"%(args.formula," ".join(scoretags[1:])))
			awkf=args.formula
			col_str=" ".join(used_tags)
			ct=1
			for tag in used_tags:
				awkf=awkf.replace(tag,'$%d'%ct)
				ct+=1
			cmd="silent_data %(infile)s %(col_str)s description | awk '{print %(awkf)s,$NF}' | sort -n -k 1 %(REVERSE)s |" \
														"head -n %(NSTRUCT)d | awk '{print $2}'"%locals()
			pipe=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
		else:
			scoretag=''
			if scorecol_defined:
        scoretag = scoretags[ abs(SCORECOL) ]
        print "use "+scoretag

			if scorecol_name_defined:
        assert( scoretags.count( scorecol_name ))
        SCORECOL = scoretags.index( scorecol_name )
        scoretagdd = scoretags[ abs(SCORECOL) ]
        scoretag = scorecol_name

			cmd="silent_data "+infile+" %(scorecol_name)s description | sort -n -k 1 %(REVERSE)s | head -n %(NSTRUCT)d | awk '{print $2}'"%locals()
			pipe=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)

		tags=[ line.strip('\n') for line in pipe.stdout ]
		count = 1
		fid = open( infile )
		line = fid.readline()

    writeout = 0
    head_score = ''
    head_sequence = ''
    while line:
        cols = string.split(line)
        if (len(cols)>1 and cols[0]=='SCORE:'):
            if tags.count(cols[-1]) > 0:
                writeout = 1
                if ( len(head_sequence) > 1 ):
                    print head_sequence
                    head_sequence=''
                if ( len(head_score) > 1 ):
                    print head_score
                    head_score=''
            else:
                writeout = 0
        if cols[1]=='score':
            writeout = 0
            head_score = line[:-1]
        if cols[0]=='SEQUENCE:':
            writeout = 0
            head_sequence = line[:-1]
        if writeout:
            print line[:-1]
        line = fid.readline()



