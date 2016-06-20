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
from os.path import exists
from os.path import basename
import argparse
from ExampleArgumentParser import ExampleArgumentParser
import sys
from amino_acids import longer_names
### toolbox library
#import library

parser = ExampleArgumentParser(prog=basename(__file__), description="obtain fasta sequences from provided pdb-files",
examples=['%(prog)s target1.pdb target2.pdb > all.fasta',
          '%(prog)s target1.pdb target2.pdb -o all.fasta'])

parser.add_argument("infiles", nargs='*',  help="pdb files");
#parser.add_argument("-nochain", dest='removechain', action='store_true')
#parser.add_argument("-chain", dest='removechain', action='store_false')
parser.add_argument('-o', metavar='all.fasta', help='output file name', default=None)
args = parser.parse_args()
if not args.infiles or len(args.infiles)==0:
    print "Require at least one pdb file as input"

pdbnames=args.infiles

if args.o:
    fastaid=open(args.o,'w')
else:
    fastaid = sys.stdout

for pdbname in pdbnames:
#    if (pdbname[-4:] != '.pdb'):
#        pdbname += '.pdb'

    outfile = pdbname

 #   removechain = args.removechain
    removechain=0
    netpdbname = pdbname
    assert( exists(netpdbname))
    #print 'Reading ... '+netpdbname

    lines = open(netpdbname,'r').readlines()

    #outid = open( outfile, 'w')
    #print 'Writing ... '+pdbname

    #fastafile = pdbname+'.fasta'
    #fastaid = open( fastafile, 'w')
    #print 'Writing ... '+fastafile
    fastaid.write('>'+basename(pdbname)+'\n');

    oldresnum = '   '
    count = 0;
    for line in lines:
        if (len(line)>20): # and (chainid == line[21]):
            line_edit = line
            if line[0:3] == 'TER':
                break
            elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit)>75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]

            if line_edit[0:4] == 'ATOM':
                resnum = line_edit[23:26]
                if not resnum == oldresnum:
                    count = count + 1
                    longname = line_edit[17:20]
                    if longer_names.has_key(longname):
                        fastaid.write( longer_names[longname] );
                    else:
                        fastaid.write( 'X')
                oldresnum = resnum

                newnum = '%3d' % count
                line_edit = line_edit[0:23] + newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+' '+line_edit[22:]

                #outid.write(line_edit)
    fastaid.write('\n')


    #outid.close()
    #fastaid.close()
