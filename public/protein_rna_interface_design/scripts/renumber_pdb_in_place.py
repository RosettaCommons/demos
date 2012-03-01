#!/usr/bin/python

import string
from sys import argv,stderr,stdout
from os import popen,system
from os.path import exists
from SWA_amino_acids import longer_names

assert( len(argv)>1)
pdbnames = argv[1:]


for pdbname in pdbnames:
	lines = open(pdbname,'r').readlines()

	oldresnum = '   '
	count = 0;

	outid  = open( 'temp.txt','w')

	atomnum  = 0
	for line in lines:
		line_edit = line
		if line[0:3] == 'TER': continue

		if line_edit[0:4] == 'ATOM' or line_edit[0:6] == 'HETATM':

			if not (line[16]==' ' or line[16]=='A'): continue

			atomnum += 1

			resnum = line_edit[23:26]
			if(resnum != oldresnum):
				count = count + 1
			oldresnum = resnum

			newnum = '%4d' % count
			line_edit = '%s%5d%s%s%s' % (line_edit[0:6],atomnum,line[11:22], newnum, line_edit[26:] )

			outid.write(line_edit)

	outid.close()

	system( 'mv temp.txt '+pdbname )
