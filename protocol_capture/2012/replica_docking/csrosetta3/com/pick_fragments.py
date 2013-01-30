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

from os import path
import shutil

import argparse
from ExampleArgumentParser import ExampleArgumentParser
import sys
import subprocess
import os
### toolbox library
import library
import traceback
import textwrap
import fasta
#############################

parser = ExampleArgumentParser(prog=path.basename(__file__),
description="Pick fragments based on chemical shifts and sequence information. This application is a wrapper that integrates several steps required to obtain fragments. First, BLAST is used to obtain a sequence profile from multiple sequence alignment. This step also provides the names of homologous proteins that can be excluded from fragment picking in benchmark mode using the -nohom flag. Second, TALOS+ is executed to obtain secondary structure predictions based on the chemical shift data. Finally, the ROSETTA application fragment_picker is started to assemble the fragment libraries using the VALL. Obviously, this wrapper has a great many dependencies. These dependencies are configured in csrosetta3/frag_picker/setup_paths.pl. Run install.py if dependencies have changed or to update the BLAST sequence database.",
examples=[('%(prog)s -cs 2jrm_trim.tab -nohom','pick fragments from non-homologous proteins (for benchmarking)'),
					('%(prog)s -cs exp.tab','pick fragments for chemical shifts in exp.tab')])
#         '%(prog)s -s 5 -e 100 full.pdb > trim.pdb'])
parser.add_argument("-cs", help="chemical shifts in TALOS format", metavar="cs.tab", default=None, required='true' )
parser.add_argument("-sizes", help="which sizes of fragments shall be build", type=int, nargs='*', default=[3, 9] )
parser.add_argument("-nfrags", help="how many frags per size-class and sequence position to collect", type=int, default=200 )
parser.add_argument("-hom", help="pick fragments also from homologous proteins [default]", action='store_true', default=True, dest='hom' )
parser.add_argument("-nohom", help="do not pick fragments from homologous proteins", action='store_false',dest='hom')
parser.add_argument("-outlier", help="report chemical shift outliers that are x*limit (no effect on fragments)", default=1.5)
parser.add_argument("-trim", help='trim the sequence according to TALOS+ output and pred2rigid', action='store_true', default=False)
parser.add_argument("-fasta", help="a target sequence", default=None )
parser.add_argument("-nocheck", dest='check', help="don't run TALOS to check for chemical shift offsets or trimmin", action='store_false', default=True)
library.add_standard_args( parser )
args = parser.parse_args()

#output:
verbose=1
library.hello(__file__)

def check_offsets(target,cs_file):
	talos_dir='%(target)s.fasta.talos'%locals()
	print 'Run TALOS+...'
	pipe=subprocess.Popen('mkdir -p %(talos_dir)s; cd %(talos_dir)s; talos+ -in ../%(cs_file)s 2>/dev/null'%locals(), shell=True, stdout=subprocess.PIPE)
	outlier_line="Checking for Chemical Shift Outliers ..."
	offset_str=""
	for line in pipe.stdout:
		line=line.strip('\n')
		tags=line.split()
		if len(tags)>4 and tags[3]=='Secondary' and abs(float(tags[5])/float(tags[7]))>args.outlier:
			if outlier_line:
				print outlier_line
				outlier_line=None
			print line
		if 'Checking Chemical Shift Referencing' in line:
			atom=None
			offset=None
			for line in pipe.stdout:
				if 'ANN prediction' in line: break
				line=line.strip('\n')
				print line
				tags2=line.split()
				if 'Referencing' in line:
					atom=tags2[4]
					if len(tags2)>5:
						offset=float(tags2[5])
				elif atom and len(tags2)>=1:
					offset=float(tags2[0])
				if atom and offset:
					if atom=='CA/CB:':
						offset_str+='CA: %(offset)f CB: %(offset)f '%locals()
					elif atom=='HA:':
						offset_str+='HA: %(offset)f '%locals()
					elif atom=="C':":
						offset_str+='C: %(offset)f '%locals()
					elif atom=="N:" or atom=="HN:" or atom=="H:":
						offset_str+='%(atom)s: %(offset)f '%locals()
					atom=None
					offset=None

	if len(offset_str):
		original_shifts=cs_file
		cs_file='%s.corrected.tab'%(path.splitext(cs_file)[0])
		print 'Re-referencing chemical shifts: the re-referenced shifts are stored to %(cs_file)s...'%locals()
		subprocess.call('correctCSoffset %(original_shifts)s -offsets "%(offset_str)s" > %(cs_file)s'%locals(), shell=True)
		print 'Rerun TALOS+ with %(cs_file)s...'%locals()
		pipe=subprocess.Popen('mkdir -p %(talos_dir)s; cd %(talos_dir)s; talos+ -in ../%(cs_file)s 2>/dev/null'%locals(), shell=True, stdout=subprocess.PIPE)
		for line in pipe.stdout:
			if 'Offset' in line:
				print 'WARNING still offsets in chemical shifts', line
	return cs_file,talos_dir

try:
	try:
		cs_root=os.environ['csrosettaDir']
	except KeyError:
		print "ERROR: Cannot find environment variable csrosettaDir."
		print "You probably forgot to run 'source <yourpath>/csrosetta3/com/init'."
		print "It is recommended to put this line into your .bashrc or .cshrc file"
		exit()


	#the perl script called from this wrapper will otherwise produce new files in a remote directory
	# this is clearly unexpected behavior hence we copy file to local dir.
	cs_file=args.cs
	if '/' in args.cs:
		if path.exists( path.basename(args.cs) ):
			raise library.RunException("Cannot overwrite %s!.\n You specified %s as the chemical shift file, but when attempting to copy this file to the local directory a file with the same name was found. Delete %s to proceed."%(path.basename(args.cs), args.cs, path.basename(args.cs)))
		print "copying chemical shifts to local directory..."
		shutil.copy(args.cs, path.basename(args.cs) )
		cs_file=path.basename(args.cs)

	target=cs_file.split('.')[0]

	frag9_file='%(target)s.frags9.dat'%locals()
	frag3_file='%(target)s.frags3.dat'%locals()
	if path.exists( frag9_file+'.gz' ) and path.exists( frag3_file+'.gz' ):
		print 'fragments exist already - remove %(frag3_file)s.gz and %(frag9_file)s.gz before running'%locals()
		exit()
	else:
		#these files will ause the fragment_pl script to skip... remove them, if the final output files are missing
		remove=['frags.fsc.score.200.9mers','frags.fsc.score.200.3mers',frag3_file+'.gz',frag9_file+'.gz']
		for r in remove:
			if path.exists(r):
				os.remove(r)


	size_flags=" ".join([ "-frag_sizes %d"%x for x in args.sizes])
	n_flag=""
	if args.nfrags:
		n_flag="-n_frags %d"%args.nfrags

	if args.hom:
		hom_flag='-hom'
	else:
		hom_flag='-nohom'


	if args.check:
		cs_file,talos_dir=check_offsets(target,cs_file)

		if args.trim:
			original_shifts=cs_file
			cs_file='%s.autotrim.tab'%(path.splitext(cs_file)[0])
			subprocess.call('cd %(talos_dir)s; pred2rigid pred.tab > pred.rigid; cd ..; renumber_talos %(original_shifts)s %(cs_file)s -rigid %(talos_dir)s/pred.rigid'%locals(), shell=True, cwd=os.getcwd())
			print 'Rerun TALOS+ with %(cs_file)s...'%locals()
			pipe=subprocess.Popen('mkdir -p %(talos_dir)s; cd %(talos_dir)s; talos+ -in ../%(cs_file)s 2>/dev/null'%locals(), shell=True, stdout=subprocess.PIPE)

	from cs import TalosCSFile
	tab=TalosCSFile()
	tab.read_file( cs_file )
	sequence=tab.sequence
	fasta_file=target+'.fasta'

	if path.exists( fasta_file ):
		target_fasta=fasta.read_fasta(fasta_file)
		if not '-' in sequence and target_fasta!=sequence:
			print "inconsistent fasta sequence: between chemical shifts %(cs_file)s and fasta file %(fasta_file)s"%locals()
			print "will overwrite fasta file, and create backup of original fasta file .bak"
			shutil.copy(fasta_file,fasta_file+".bak")
			fasta.write_fasta(fasta_file, sequence, target )

	fasta_file=target+'.fasta'
	if '-' in sequence and args.fasta:
		target_fasta=fasta.read_fasta(args.fasta)
		if not '-' in target_fasta:
			print 'missing residues in the sequence information of the talos-file, will use sequence from fasta file %s'%fasta_file
			start=fasta.find_fasta_offset(target_fasta,sequence)
			end=start+len(sequence)-1;
			if start!=0:
				print 'WARNING: using fasta sequence to fill gaps, but will use start and end from the talos file. Use renumber_talos to fix this before running frag_picker'
			if start>=0:
				sequence=target_fasta[start:end]
			else:
				raise library.MissingInput('provided sequence is trimmed to much at N-terminus')
			tab.sequence=sequence
			print sequence
			tab.write_file( cs_file )
			print cs_file
		else:
			exit('require full sequence information, provide fasta file with -fasta')


	fasta.write_fasta(fasta_file, sequence, target )

	subprocess.call(('perl %(cs_root)s/frag_picker/pick_fragments.pl -talos_fn %(cs_file)s -fasta %(fasta_file)s %(hom_flag)s %(size_flags)s %(n_flag)s'+
									'&& mv frags.score.200.9mers %(frag9_file)s && mv frags.score.200.3mers %(frag3_file)s && gzip %(target)s.frags?.dat')%locals(), shell=True, cwd=os.getcwd())

except library.LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)

