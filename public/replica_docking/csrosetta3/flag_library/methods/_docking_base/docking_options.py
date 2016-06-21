##-*- mode:python;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments

import string
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen,mkdir,makedirs
from os import dup2,path
from os.path import exists
from os.path import *
from os.path import basename
import argparse
import sys
import shutil
### toolbox library
from library import Tracer
import library
import automatic_setup

#from _presetup import TargetLib
import traceback

# definition of options for the your method

if 'group' in locals():
	group.add_argument("-native", help="native complex" );
	group.add_argument("-prepack_native", help="native complex prepacked with unbound structures' sidechains" );
	group.add_argument("-pdb", help="start conformation randomized on the native" );
	group.add_argument("-partners", help="docking partners" );
	group.add_argument("-cst_file", help="cst_file" );
#	group.add_argument("-disulf", help="disulfide bonds pairs");

if 'run_group' in locals():
	run_group.add_argument("-score", help="score weights, options: interchain_cen, score3, docking, score12..." );
	run_group.add_argument("-extra_score", choices=['docking_interface_score'], help="" );
	run_group.add_argument("-min_score", help="hack on the low resolution score function")
	run_group.add_argument("-nstruct", help="of each batch if batches used, default=1", default=1 );
	run_group.add_argument("-start", help="path to decoy-sets or pdb file for refinement (absolut)" );
	run_group.add_argument("-cst_weight", help="default=5", default=5)



tr = Tracer( "Dock" )
# method-based code
# DockMethod is the base-clase for both rosetta_dock and replica_dock

class DockMethod(automatic_setup.BasicMethod):
	def __init__(self,name,path):
		automatic_setup.BasicMethod.__init__(self,name,path)
		self.non_file_options.append('score')
		self.non_file_options.append('extra_score')
		self.non_file_options.append('min_score')
		self.non_file_options.append('nstruct')
		self.non_file_options.append('start')
		self.non_file_options.append('cst_weight')

		self.option2dir['native']='native'
		self.option2dir['pre_native']='start'
		self.option2dir['pdb']='start'
		self.option2dir['partners']='partners'
		self.option2dir['cst_file']='cst'

	def make_target_flags(self, run, setup, filename, flags, subs ):
		tr.out("make target flags...")
		args=self.get_args()

		# parameters in common for different docking protocol

# 		if args.partners:
# 			run.add_subst('CM_PARTNERS', setup.cm_path(args.partners))
# 			flags.write("@$CM_PARTNERS\n")

		# score
# 		if args.score:
# 			flags.write("-score:weights %s\n"%args.score)

		# extra_score
		if args.extra_score:
			flags.write("-score:docking_interface_score")

		# min_score
		if args.extra_score:
			flags.write("-score:min_score_score %f\n"%float(args.min_score))

# 		# cst   # transfered to rosetta_dock/options.py
# 		if args.cst_file:
# 			flags.write("-cst_file %s\n"%setup.cm_path(args.cst_file))
# 			flags.write("-cst_weight %s\n"%args.cst_weight)

		# writing everything left into flags
		flags.write("-nstruct %d\n"%int(args.nstruct))

		# native and docking partners(for which only multi-chain protein is neccessary)
		if args.native:
			run.add_subst('CM_NATIVE_PDB', setup.cm_path(args.native))
			flags.write("-in:file:native $CM_NATIVE_PDB\n")

