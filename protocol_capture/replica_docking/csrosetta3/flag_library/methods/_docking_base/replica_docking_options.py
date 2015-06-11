##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:3 -*-'
## make mammoth structure alignments

import string
from glob import glob
from os import dup2,path
from os.path import exists
from os.path import basename
import argparse
import sys
import shutil

### toolbox library
from library import Tracer
import library
from automatic_setup import BasicMethod
import traceback

sub_method_code = flag_lib+"/methods/_docking_base/docking_options.py"
if exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()

# definition of options for the method
# if 'group' in locals():
# 	 group.add_argument("-partners", help="partner-string to specify which chains are docked", default="A_B" );

if 'run_group' in locals():
	 run_group.add_argument("-protocol", choices=['rep_high','rep_cen','rep_high_rigid','rep_high_frozen'], help="choose the mode", default="rep_cen");
	 run_group.add_argument("-xml", choices=['u_inter', 'uniform','gaussian'], help="u_inter: unbiased rotation sampling with interpolated step_size;\nuniform: unbiased rotation sampling with constant step_size;\ngaussian: original gaussian_move, biased rotation sampling");
	 run_group.add_argument("-n_replica", help="3: normal replica dock centroid;\n7: replica dock centroid with min_score", default=3);

tr = Tracer( "replica_docking_base_method" )

# method-based code
class ReplicaDockingBaseMethod(DockMethod):
	 def __init__(self,name,path):
			DockMethod.__init__(self,name,path)
			self.non_file_options.append('protocol')
			self.non_file_options.append('xml')

	 def make_target_flags(self, run, setup, filename, flags, subs ): # for the non-fixed options
			DockMethod.make_target_flags( self, run, setup, filename, flags, subs )
			tr.out("make target flags...")
			args=self.get_args()

			if args.partners:
				 run.add_subst('CM_PARTNERS', setup.cm_path(args.partners))
				 flags.write("@$CM_PARTNERS\n")

			# start file
			if args.start: # if give a pdb file
				 if isfile(args.start):
						fname,fext = splitext(args.start)
						if fext==".pdb":
							 run.add_subst('CM_INPUT_PDB', args.start)
				 print "only accept .pdb file for -start in replica dock, will use default file P_*.pdb from the target library instead of your unacceptable %s\n"%args.start
		# no given pdb file, then use default pdb files from the target library
			else:
				 if not args.pdb:
						raise library.MissingInput("P.pdb does not exist, you need a start pdb for replica dock")
				 run.add_subst('CM_INPUT_PDB', setup.cm_path(args.pdb))
			flags.write("-in:file:s $CM_INPUT_PDB\n")

#			flags.write("-nstruct %d\n"%int(args.nstruct) )

			if args.protocol == "rep_cen":
				 if args.score == "interchain_cen":
						flags.write("-score:weights %s\n"%args.score)
				 else:
						raise library.MissInput("wrong score weights for low resolution docking")
				 if args.min_score:
						flags.write("-score:min_score_score %f\n"%float(args.min_score))
			else:
				 if not args.score or args.score=="interchain_cen":
						raise library.MissInput("you need to give score weights correctly for high res docking")
				 flags.write("-score:weights %s\n"%args.score)
				 if args.extra_score:
						flags.write("-score:docking_interface_score")

			if not args.n_replica:
				 raise library.MissingInput("you need to specify n_replica for method %s"%self.name)
			flags.write("-n_replica %d\n"%int(args.n_replica))


	 def setup_file_library( self ):
 			BasicMethod.setup_file_library( self )
			args = self.get_args()
			fl = self.file_library
			path = flag_lib+"/methods/_docking_base/"

			fl.executable  = "rosetta_scripts"

			fl.provide_file( "flags", path, "flags_docking_replica" )
			fl.add_string( "commandline", "@flags_docking_replica @$CM_FLAGFILE ")

			### the following part could always be improved!
			if args.protocol == 'rep_cen':
				 if int(args.n_replica) == 7:
						fl.provide_file( "flags", path, "hamiltonians_cen_7.txt")
						if args.xml=='u_inter':
							 fl.provide_file("flags", path, "dock_cen_inter_7.xml")
							 fl.add_string( "commandline", " -parser:protocol @@dock_cen_inter_7.xml ")
						if args.xml=='uniform':
							 fl.provide_file("flags", path, "dock_cen_7.xml")
							 fl.add_string( "commandline", " -parser:protocol @@dock_cen_7.xml ")
						if args.xml=='gaussian':
							 fl.provide_file("flags", path, "dock_cen_g_7.xml")
							 fl.add_string( "commandline", " -parser:protocol @@dock_cen_g_7.xml ")
				 if int(args.n_replica) == 3:
						fl.provide_file( "flags", path, "hamiltonians_cen.txt")
				 #				 fl.provide_file( "flags", path, "exchange_schedule_cen.txt")
						if args.xml=='u_inter':
							 fl.provide_file("flags", path, "dock_cen_inter.xml")
							 fl.add_string( "commandline", " -parser:protocol @@dock_cen_inter.xml ")
						if args.xml=='uniform':
							 fl.provide_file("flags", path, "dock_cen.xml")
							 fl.add_string( "commandline", " -parser:protocol @@dock_cen.xml ")
						if args.xml=='gaussian':
							 fl.provide_file("flags", path, "dock_cen_g.xml")
							 fl.add_string( "commandline", " -parser:protocol @@dock_cen_g.xml ")
			if args.protocol == 'rep_high':
				 fl.provide_file( "flags", path, "hamiltonians_high.txt")
				 if args.xml=='u_inter':
						fl.provide_file("flags", path, "dock_u_inter.xml")
						fl.add_string( "commandline", " -parser:protocol @@dock_inter.xml ")
				 if args.xml=='uniform':
						fl.provide_file("flags", path, "dock_u.xml")
						fl.add_string( "commandline", " -parser:protocol @@dock.xml ")
				 if args.xml=='gaussian':
						fl.provide_file("flags", path, "dock_g.xml")
						fl.add_string( "commandline", " -parser:protocol @@dock_g.xml ")
			if args.protocol == 'rep_high_rigid':
				 fl.provide_file( "flags", path, "dock_high_rigid.xml" )
				 fl.add_string( "commandline", "-parser:protocol @@dock_high_rigid.xml ")
			if args.protocol == 'rep_high_frozen':
				 fl.provide_file( "flags", path, "dock_high_frozen.xml" )
				 fl.add_string( "commandline", "-parser:protocol @@dock_high_frozen.xml ")


