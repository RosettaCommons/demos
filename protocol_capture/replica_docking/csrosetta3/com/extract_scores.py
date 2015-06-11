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

import string
from silent_lib import ReadSilentData
import sys

from os.path import basename
import argparse
from ExampleArgumentParser import ExampleArgumentParser
import traceback
import library


#avoid problem with broken pipes (e.g., head)
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

parser = ExampleArgumentParser(prog=basename(__file__), description="obtain fasta sequences from provided pdb-files",
examples=['%(prog)s silent_file.out rms score > rms_score.txt',
          '%(prog)s scores_silent.fsc rms atom_pair_constraint rdc > rms_score.txt'],
aliases=['extract_scores','silent_data'])

parser.add_argument("silentfile", help="input silent file", default=None);
parser.add_argument('columns', nargs='*', help='columns for output');
parser.add_argument('-pipe', help="read from stdin", default=False, action='store_true');

library.add_standard_args( parser )

args = parser.parse_args()
names = args.columns

try:
  if args.pipe:
    infile=sys.stdin
    names=[args.silentfile]+args.columns
  else:
    infile = open(args.silentfile,'r')
  sfd = ReadSilentData( names )
  for l in infile:
    data=sfd.read_line( l )
    if data:
      if isinstance( data, str ):
        print data
      else:
        print " ".join(data)

except library.LibException as inst:
  if args.traceback:
    print traceback.print_exc(inst )
  else:
    print traceback.print_exception(sys.exc_type, sys.exc_value, None)
