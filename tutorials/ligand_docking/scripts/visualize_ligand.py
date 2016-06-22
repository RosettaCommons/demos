#!/usr/bin/env python

# system stuff
import os
import copy
from optparse import OptionParser
import sys
import subprocess
#from sys import argv

# for importing as a plugin into PyMol
import __main__
__main__.pymol_argv = ['pymol', '-qc']
import pymol
from pymol import cmd
from pymol import stored
from pymol import selector

pymol.finish_launching()
pymol.cmd.feedback('enable', 'all', 'actions')
pymol.cmd.feedback('enable', 'all', 'results')

from array import *

### Parse the options ###
parser=OptionParser(usage='visualize_ligand.py [PDB-FILE-1] [PDB-FILE-2] ...')

from sys import argv
if len(argv) == 1 : # only the script name
	argv.append('--help')

(options, args)=parser.parse_args()

for file in args:
	pymol.cmd.delete('docked')
	pymol.cmd.load(file, 'docked')

	# make selections of ligand, binding pocket (all residues with an atom within 4A of ligand) and thei rpolar contacts
	pymol.cmd.select(name='ligand', selection='docked and hetatm')
	pymol.cmd.select(name='pocket', selection='ligand around 4 and docked')
	pymol.cmd.select(name='pocket', selection='byres pocket')
	pymol.cmd.distance('contacts', 'ligand', 'pocket', 4, 2)
	pymol.cmd.color(4,'contacts')
		
	#stylistic changes
	pymol.preset.ligand_cartoon('all')
	pymol.cmd.delete('docked_pol_conts')
	pymol.cmd.show('cartoon','all')
	pymol.cmd.hide('lines')
	pymol.cmd.bg_color('black')
	pymol.cmd.set('cartoon_fancy_helices', 1)
	pymol.cmd.show('sticks','pocket & (!(n;c,o,h|(n. n&!r. pro)))')
	pymol.cmd.show('sticks','ligand')
	pymol.cmd.hide('lines', 'hydro')
	pymol.cmd.hide('sticks', 'hydro')
	
	#save session
	filename=os.path.splitext(file)[0]
	outfile=filename+'.pse'
	pymol.cmd.save(outfile)	










