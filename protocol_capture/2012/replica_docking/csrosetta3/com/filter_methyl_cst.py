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
## -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-

#application specific headers

#default headers
import argparse
from ExampleArgumentParser import ExampleArgumentParser
from os.path import basename
import traceback, sys

#toolbox headers
import cs
import library
import amino_acids

parser = ExampleArgumentParser(prog=basename(__file__),
                                 fromfile_prefix_chars='@',
                                 description="Remove restraints that use side-chain atoms of unallowed amino-acid types. "+
															 "This application provides a quick way to simulate, e.g.,an ILV-labelled restraints file from a "+
															 "fully labelled restraints file.",
															 examples=[('%(prog)s -pdb ref.pdb -allowed ILV all.cst ilv_only.cst',
																					'copy all restraints from all.cst into ilv_only.cst, which do not contain side-chain atoms from other residues'),
																				 ('%(prog)s -fasta ref.fasta -allowed ILVAMT all.cst ilv_only.cst',
																					'copy all restraints from all.cst into ilv_only.cst, which do not contain side-chain atoms from other residues')],
                          )

parser.add_argument("-pdb", help="pdb file to get sequence");
parser.add_argument("-fasta", help="fasta file");
parser.add_argument("-allowed", help="string of allowed amino-acids such as ILV", default="ILV");
parser.add_argument("cstin", help="rosetta restraint file");
parser.add_argument("cstout", help="rosetta restraint file",nargs="?");
parser.add_argument("-resfile", help="output the residue-numbers of allowed sidechains");

args = parser.parse_args()

#output:
if args.cstout:
	verbose=1
	cstout=open(args.cstout,'w');
else:
	cstout=sys.stdout
	verbose=0
####### program start
if verbose:
	library.hello( __file__ )

fasta=0
if args.fasta:
	fasta=library.read_fasta(args.fasta)
elif args.pdb:
	fasta=library.pdb2fasta(args.pdb);
else:
	print
	print "ERROR: either fasta or pbd file are required to determine residue numbers for filtering...\n"
	exit()

resi=1;
sidechains=[];
for aa in fasta:
    if aa in args.allowed:
        sidechains.append( resi );
    resi=resi+1;

if verbose: print sidechains;
if args.resfile:
    resfile=open(args.resfile,'w');
    for i in sidechains:
        resfile.write("%d\n"%i);

sidechains.append( 0 ) #0 is used for backbone hydrogen atoms --- allways allowed
fasta.insert(0,"X");
lines=open( args.cstin, 'r').readlines();
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
    bad = 2
    if (l[1] == 'H'):
        resi1 = 0
    if ( l[3] == 'H'):
        resi2 = 0
#    print res_lines
    if resi1 in sidechains and resi2 in sidechains:
        cstout.write(line);
    else:
        cstout.write("# %c-%c IS NOT %s restraint --- removed: %s"%( fasta[resi1],fasta[resi2],args.allowed, line))

