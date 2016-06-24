#!/usr/bin/env python

from optparse import OptionParser
import sys
import subprocess
from sys import argv

from array import *
 
# system stuff
import os
import copy
 
# pretty printing
import pprint
 
# for importing as a plugin into PyMol
 
# using numpy for linear algebra
import numpy

### Parse the options ###
parser=OptionParser(usage='calculate_rmsd.py [options] <list of pdbs>', description="\
This script calculates RMS to native for each provided model.\
It does so for 1) the whole complex, 2) all HETATMs or specified chains,\
3) '2' + side-chains within specified angstroms\n\
")

from sys import argv
if len(argv) == 1 : # only the script name
	argv.append('--help')

### Add options
parser.add_option("-n", "--native", dest="native", type='string',
		help="path to the native PDB structure")
parser.add_option("-c", "--chains", dest="chains", type='string',
		help="comma separated list of 1-char chains as found in PDB")
parser.add_option("-a", "--angstroms", dest="angstroms", type='string', default="4",
		help="protein residues within this many angstroms will be included in a separate RMS calculation")
parser.add_option("-o", "--output", dest="output", type='string', default="out.txt",
		help="path to output file")
parser.add_option("--include_hydrogens", dest="include_hydrogens", default="False", action="store_true",
		help="provide this flag to include Hs in the calculations")
parser.add_option("-s", "--scores_to_keep", dest="scores_to_keep", type='string', #default="total_score,interface_delta_X,ligand_rms_no_super_X",
		help="comma separated list of score terms from the PDB to include in the table")
parser.add_option("--symmetrical_chains", dest="symmetrical_chains", type='string',
		help="comma separated list of protein chains that are super_imposeable. This script will report the minimum RMS from aligning each pairing of superimposeable chains")
parser.add_option("--no_align", dest="do_align", default="True", action="store_false", #default does align before RMSD
		help="provide this flag to calculate RMSD without doing alignment")

(options, args)=parser.parse_args()

import __main__
__main__.pymol_argv = ['pymol', '-qc']
import pymol
from pymol import cmd
from pymol import stored
from pymol import selector

pymol.finish_launching()
pymol.cmd.feedback('disable', 'all', 'actions')
pymol.cmd.feedback('disable', 'all', 'results')

import gzip
def find_scores(file):
	lines=None
	if file.endswith('.gz'):
		lines=gzip.open(file).readlines()
	else:
		lines=open(file).readlines()
	scores=[]
	for score_term in options.scores_to_keep.split(','):
		for line in reversed(lines):
			if line.startswith('#'):# assumes all score terms are at the end of the file, no comments after them
				sys.exit(score_term+" was not found in "+file)
			(this_term, score)=line.split()[:2]
			if this_term==score_term:
				scores.append(score)
				break
	return scores

def align_for_min_rms(file):
	print "Doing aligning"
	if options.symmetrical_chains==None:
		pymol.cmd.align('docked and not hetatm', 'native and not hetatm')
		return
	#else
	best_rms=100
	best_alignement=(None,None)
	symmetrical_chains=options.symmetrical_chains.split(',')
	for i in range(len(symmetrical_chains)-1):
		for j in range(i+1, len(symmetrical_chains)):
			pymol.cmd.align('docked and not hetatm and chain '+symmetrical_chains[i], 'native and not hetatm and chain '+symmetrical_chains[j])
			this_rms= calculate_ligand_rms()
			if this_rms < best_rms:
				best_rms = this_rms
				best_alignment=(i,j)
	pymol.cmd.align('docked and chain '+symmetrical_chains[best_alignment[0]], 'native and chain '+symmetrical_chains[best_alignment[1]])

def calculate_ligand_rms():
		pymol.cmd.rms_cur('docked_lig', 'native_lig')
def optAlign( sel1, sel2 ):
	"""
	optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
	Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
	Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA
 
	Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
	PyMol based upon your selections.
 
	By default, this program will optimally align the ALPHA CARBONS of the selections provided.
	To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
 
	@param sel1: First PyMol selection with N-atoms
	@param sel2: Second PyMol selection with N-atoms
	"""
	cmd.reset()
 
	# make the lists for holding coordinates
	# partial lists
	stored.sel1 = []
	stored.sel2 = []
	# full lists
	stored.mol1 = []
	stored.mol2 = []
 
	# Get the selected coordinates.  We
	# align these coords.
	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
 
	# get molecule name
	mol1 = cmd.identify(sel1,1)[0][0]
	mol2 = cmd.identify(sel2,1)[0][0]
 
	# Get all molecule coords.  We do this because
	# we have to rotate the whole molcule, not just
	# the aligned selection
	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
 
	# check for consistency
	assert len(stored.sel1) == len(stored.sel2)
	L = len(stored.sel1)
	assert L > 0
 
	# must alway center the two proteins to avoid
	# affine transformations.  Center the two proteins
	# to their selections.
	COM1 = numpy.sum(stored.sel1,axis=0) / float(L)
	COM2 = numpy.sum(stored.sel2,axis=0) / float(L)
	stored.sel1 -= COM1
	stored.sel2 -= COM2
 
	# Initial residual, see Kabsch.
	E0 = numpy.sum( numpy.sum(stored.sel1 * stored.sel1,axis=0),axis=0) + numpy.sum( numpy.sum(stored.sel2 * stored.sel2,axis=0),axis=0)
 
	#
	# This beautiful step provides the answer.  V and Wt are the orthonormal
	# bases that when multiplied by each other give us the rotation matrix, U.
	# S, (Sigma, from SVD) provides us with the error!  Isn't SVD great!
	V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(stored.sel2), stored.sel1))
 
	# we already have our solution, in the results from SVD.
	# we just need to check for reflections and then produce
	# the rotation.  V and Wt are orthonormal, so their det's
	# are +/-1.
	reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
 
	if reflect == -1.0:
		S[-1] = -S[-1]
		V[:,-1] = -V[:,-1]
 
	RMSD = E0 - (2.0 * sum(S))
	RMSD = numpy.sqrt(abs(RMSD / L))
 
	#U is simply V*Wt
	U = numpy.dot(V, Wt)
 
	# rotate and translate the molecule
	stored.sel2 = numpy.dot((stored.mol2 - COM2), U)
	stored.sel2 = stored.sel2.tolist()
	# center the molecule
	stored.sel1 = stored.mol1 - COM1
	stored.sel1 = stored.sel1.tolist()
 
	# let PyMol know about the changes to the coordinates
	cmd.alter_state(1,mol1,"(x,y,z)=stored.sel1.pop(0)")
	cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
 
	# we already have our solution, in the results from SVD.
	# we just need to check for reflections and then produce
	# the rotation.  V and Wt are orthonormal, so their det's
	# are +/-1.
	reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
 
	if reflect == -1.0:
		S[-1] = -S[-1]
		V[:,-1] = -V[:,-1]
 
	RMSD = E0 - (2.0 * sum(S))
	RMSD = numpy.sqrt(abs(RMSD / L))
 
	#U is simply V*Wt
	U = numpy.dot(V, Wt)
 
	# rotate and translate the molecule
	stored.sel2 = numpy.dot((stored.mol2 - COM2), U)
	stored.sel2 = stored.sel2.tolist()
	# center the molecule
	stored.sel1 = stored.mol1 - COM1
	stored.sel1 = stored.sel1.tolist()
 
	# let PyMol know about the changes to the coordinates
	cmd.alter_state(1,mol1,"(x,y,z)=stored.sel1.pop(0)")
	cmd.alter_state(1,mol2,"(x,y,z)=stored.sel2.pop(0)")
 
cmd.extend("optAlign", optAlign)

def load_file(file, name):
	print "file, name: "+file+' '+name
	pymol.cmd.load(file, name)
	pymol.cmd.alter(selection='resn TP3', expression="resn='WAT'")# call all waters by the same name

def get_hydro_string():
	if options.include_hydrogens:
		return 'and not hydro'
	return ''

def make_selections(file):
	hydro=get_hydro_string()
	load_file(options.native, 'native')
	load_file(file, 'docked')
	# make selections
	if options.chains==None:
		pymol.cmd.select(name='docked_lig', selection='docked and hetatm '+hydro)
		pymol.cmd.select(name='native_lig', selection='native and hetatm '+hydro)
	else:
		chains= options.chains.replace(',','+')
		pymol.cmd.select(name='docked_lig', selection='docked and chain '+chains+' '+hydro)
		pymol.cmd.select(name='native_lig', selection='native and chain '+chains+' '+hydro)
	pymol.cmd.select(name='docked_lig_sc', selection='(docked within '+options.angstroms+' of docked_lig) '+hydro)
	pymol.cmd.select(name='native_lig_sc', selection='(native within '+options.angstroms+' of native_lig) '+hydro)
	#raw_input()

def calculate_rms(file):
	make_selections(file)
	if options.do_align:
		align_for_min_rms(file)
	all_rms= pymol.cmd.rms_cur('docked', 'native')
	ligand_rmss=[]
	lig_and_sc_rms= pymol.cmd.rms_cur('docked_lig_sc', 'native_lig_sc')
	if options.chains==None:
		ligand_rmss.append(  pymol.cmd.rms_cur('docked_lig', 'native_lig') )
	else:
		chains= options.chains.replace(',','+')
		hydro=get_hydro_string()
		for chain in chains.split('+'):
			ligand_rms= pymol.cmd.rms_cur('docked and chain '+chain+' '+hydro, 'native and chain '+chain+' '+hydro)
			#raw_input()
			ligand_rmss.append(ligand_rms)
			#Align ligand then do RMSD
			optAlign('docked_lig', 'native_lig')
			align_ligand_rms= pymol.cmd.rms_cur('docked_lig', 'native_lig')
			ligand_rmss.append(align_ligand_rms)
		if len(chains)>1:
			ligand_rmss.append( pymol.cmd.rms_cur('docked_lig', 'native_lig') )
	rms_line=file+' '+str(round(all_rms, 2))+' '
	for rms in ligand_rmss:
		rms_line += str(round(rms, 2))+' '
	rms_line += str(round(lig_and_sc_rms, 2))
	pymol.cmd.delete('all')
	return rms_line

def get_header():
	header='File all_rms '
	if options.chains==None:
		header+='hetatm_rms '
	else:
		chains= options.chains.split(',')
		for chain in chains:
			header+=chain+'_rms '
		if(len(chains)>1):
			header+='all_chain_rms '
	header+='super_rms '	
	header+='lig_sc_rms '
	if options.scores_to_keep is not None:
		terms=options.scores_to_keep.split(',')
		for term in terms:
			header+=term+' '
	return header

output=open(options.output, 'a')
output.write( get_header()+'\n' )
for file in args:
	output.write(calculate_rms(file)+' ')
	if options.scores_to_keep is not None:
		for score in find_scores(file):
			output.write(score+' ')
	output.write('\n')
output.close()


subprocess.call(["capture_command.sh"] + sys.argv)
