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
if 'group' in locals():
	group.add_argument("-disulf", help="disulfide bonds pairs");
if 'run_group' in locals():
	run_group.add_argument("-protocol", choices=['standard','centroid','refine','dock_pert']);
#	run_group.add_argument("-disulf", help="disulfide bonds pairs");
	run_group.add_argument("-batches", help="total decoys = nstruct*batches if use" );
	run_group.add_argument("-suffix", help="suffix for decoys' description");
	run_group.add_argument("-out", help="output file name", default="dummy.out" );
	run_group.add_argument("-pattern",help="name pattern of the start decoys",default="decoys_0*.out");
	run_group.add_argument("-prefix", help="prefix for output filename based on input filename", default="refine")

tr = Tracer( "RosettaDock" )
# method-based code
sub_method_code = flag_lib+"/methods/_docking_base/docking_options.py"
if exists( sub_method_code ):
	exec open( sub_method_code, 'r' )
else:
	print "CANNOT FIND METHOD CODE %s"%sub_method_code
	exit()

class RosettaDockMethod(DockMethod):
	def __init__(self,name,path):
		print "RosettaDockMethod: ", name
		DockMethod.__init__(self,name,path)
		self.non_file_options.append('protocol')
		self.option2dir['disulf']='disulf'

	def make_target_flags(self, run, setup, filename, flags, subs ):
		DockMethod.make_target_flags( self, run, setup, filename, flags, subs )
		args=self.get_args()

		# file name of output
		if args.out:
# 			if not args.batches:
# 				if not args.protocol=="refine":
# 					flags.write("-out:file:silent %s\n"%args.out)
# 				else:
# 					if not args.start or not isdir(args.start):
# 						flags.write("-out:file:silent %s\n"%args.out)
			flags.write("-out:file:silent %s\n"%args.out)
		# score
		if args.score:
			flags.write("-score:weights %s\n"%args.score)

		# special flags for refine
		if not args.protocol=="centroid":
			#flags.write("-docking:recover_sidechains $CM_NATIVE_PDB\n")
			flags.write("-use_input_sc\n")
			flags.write("-unboundrot $CM_NATIVE_PDB\n")
		if args.protocol=="refine":   # I have a suspicious that these disulf will work automatically in the standard/dock_pert mode, I mean it probably only breaks when doing the two low/high res phase separately
			if args.disulf:
				flags.write("\n-detect_disulf true\n")
				flags.write("-rebuild_disulf true\n")
				flags.write("-fix_disulf %s\n\n"%setup.cm_path(args.disulf))
		else:
			if args.cst_file:   # we don't need cst_file for refinement
				flags.write("-cst_file %s\n"%setup.cm_path(args.cst_file))
				flags.write("-cst_weight %s\n"%args.cst_weight)

		# input or call it start file
		if not args.start: # not input, then use the default files as start
			if args.protocol == "standard" or args.protocol=="centroid":
				run.add_subst('CM_INPUT_PDB', setup.cm_path(args.pdb))
			if args.protocol == "dock_pert":
				run.add_subst('CM_INPUT_PDB', setup.cm_path(args.prepack_native))
			if args.protocol=="refine":
				run.add_subst('CM_INPUT_PDB', setup.cm_path(args.native))
			flags.write("-in:file:s $CM_INPUT_PDB\n")
			if args.partners:
				run.add_subst('CM_PARTNERS', setup.cm_path(args.partners))
				flags.write("@$CM_PARTNERS\n")
		else:
			if isfile(args.start):
				fname,fext = splitext(args.start)
				if fext=='.pdb':
					run.add_subst('CM_INPUT_PDB', args.start)
					flags.write("-in:file:s $CM_INPUT_PDB\n")
					if args.partners:
						run.add_subst('CM_PARTNERS', setup.cm_path(args.partners))
						flags.write("@$CM_PARTNERS\n")
				elif fext=='.out':
					flags.write("-in:file:silent %s\n"%args.start)
					flags.write("-out:file:silent refine_%s\n"%(basename(args.start)))
				else:
					raise library.MissingInput("file input expect .pdb or .out format")
			else:
				if not isdir(args.start):
					raise library.MissingInput("non-pdb input expected to be an absolut path to decoys-sets")
				if not args.protocol == "refine":
					raise library.MissingInput("only refine accept decoys as input right now")
				path_to_input=args.start.replace('TARGET', setup.target.name)
				assert path_to_input[0]=='/', 'absolute path required for option args.start when it is not pdb file'
				print path_to_input
				import glob
				input_decoys=glob.glob(path_to_input+'/'+args.pattern)
				print input_decoys
				batches=[]
				for i, file in enumerate(input_decoys):
					flag_file=setup.create_file("flags_batch%04d"%(i+1))
					if args.suffix:
						open(setup.abspath(flag_file),'w').write('-in:file:silent %s\n-out:file:silent %s_%s\n-out::suffix %s\n'%(file,args.prefix,basename(file),args.suffix))
					else:
						open(setup.abspath(flag_file),'w').write('-in:file:silent %s\n-out:file:silent %s_%s\n'%(file,args.prefix,basename(file)))
					batches.append(setup.abspath(flag_file))
				flags.write("\n-run:batches "+" ".join(batches)+"\n")

		# parameters which are common for bound/unbound docking
		# run batches
		if args.batches and not args.start:
			batches=[]
			for i in range(0, int(args.batches)):
				flag_file=setup.create_file("flags_batch%04d"%(i+1))
				open(setup.abspath(flag_file),'w').write("-out:file:silent decoys_%04d.out\n-out::suffix _%04d\n"%( (i+1),(i+1)) )
				batches.append(setup.abspath(flag_file))
			flags.write("\n-run:batches "+" ".join(batches)+"\n")


   def setup_file_library( self ):
		automatic_setup.BasicMethod.setup_file_library( self )
		fl = self.file_library
		args=self.get_args()
		path = flag_lib+"/methods/rosetta_dock/"

		fl.executable = "docking_protocol"

		if args.protocol == "refine":
			fl.provide_file( "flags", path, "flags_docking_refine" )
			fl.add_string("commandline", "@flags_docking_refine @$CM_FLAGFILE");
		if args.protocol == "standard":
			fl.provide_file( "flags", path, "flags_docking_standard" )
			fl.add_string( "commandline", "@flags_docking_standard @$CM_FLAGFILE");
		if args.protocol == "centroid":
			fl.provide_file( "flags", path, "flags_docking_centroid" )
			fl.add_string( "commandline", "@flags_docking_centroid @$CM_FLAGFILE");
		if args.protocol == "dock_pert":
			fl.provide_file( "flags", path, "flags_docking_pert" )
			fl.add_string( "commandline", "@flags_docking_pert @$CM_FLAGFILE");

		# fl.provide_file( "flags", path, "flags_nmr_patches" )
		# fl.override("flags", "flags_denovo", '-increase_cycles %f'%float(args.cycle_factor))
		# fl.add_string( "commandline", "-out:file:silent decoys.out @flags_denovo @$CM_FLAGFILE");

method = RosettaDockMethod(method_name, method_path)
