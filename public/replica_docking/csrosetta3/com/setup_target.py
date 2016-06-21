#!/usr/bin/env python2.7
##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
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

import argparse
import sys
import traceback

### toolbox library
import automatic_setup
from library import LibException, MissingInput
import library



#setup common-cmdline options
application = automatic_setup.SetupApplication( 'manipulate a specific setup of a target',__file__, sys.argv, examples=[
		('%(prog)s -target t000_ -method abrelax -frags frags9.dat frags3.dat -fasta t000_.fasta',
		 'prepare a basic abrelax run for target t000_ (use default target-library location $CS3_BENCH_TARGETLIB)'),
		('%(prog)s -target t000_ -method rasrec -transfer_method abrelax',
		 'copy settings from Setup (t000_, abrelax, standard) to Setup (t000_, rasrec, standard)'),
		('%(prog)s -target t000_ -method rasrec -transfer_label standard -label with_rdc -rdc med1.rdc med2.rdc',
		 'add rdc data to existing Setup (t000_, rasrec, standard) and call it Setup (t000_, rasrec, with_rdc)')]

																)

#setup specific-cmdline options
parser=application.parser
parser.add_argument("-target", help="target_dir is target_prefix/target", default="t000_" );
parser.add_argument("-transfer_method",help="start with this setup and modify according to flags", default=None )
parser.add_argument("-transfer_label", help="start with this setup and modify according to flags", default=None )
parser.add_argument("-overwrite", help="force overwriting of files in the target-database", action='store_true', default=False )
library.add_standard_args( parser )
#check usage  -- a method must be chosen
application.exit_if_no_method_selected()
args=parser.parse_args()

#main program
try:

	#setup target
	target=application.get_target( args.target )
	target.overwrite=args.overwrite
	print target

	#setup input setup  -- method and label can be chosen to be different with transfer_XXX
	input_method_name=application.method.name
	if args.transfer_method:
		input_method_name=args.transfer_method

	input_label=args.label
	if args.transfer_label:
		input_label=args.transfer_label

	input_setup = automatic_setup.Setup( target, input_method_name, input_label )
	print input_setup

	#load options into input setup
	args=application.load_setup(input_setup)

	#setup output setup
	output_setup = automatic_setup.Setup( target, application.method.name, args.label )
	if input_setup.name!=output_setup.name and not input_setup.exists():
		raise MissingInput("Setup '%s' does not exist!"%input_setup.name)
	print output_setup

	#modify setup using cmdline and store
	application.modify_setup( target, output_setup )


except LibException as inst:
	if args.traceback:
		print traceback.print_exc(inst )
	else:
		print traceback.print_exception(sys.exc_type, sys.exc_value, None)
	exit(1)


#parser.add_argument("-frags", nargs=2, help="fragment files");
#parser.add_argument("-fasta", help="fasta file");
#parser.add_argument("-target", help="target_dir is target_prefix/target", default="t000_" );
#parser.add_argument("-target_dir", help="input-files are used from/put into this path, instead of directly into the run-directory");
#parser.add_argument("-target_prefix", help="target dirs are found here" );
#parser.add_argument("-run_dir", help="setup a run-directory here");
#parser.add_argument("-quiet", help="suppress output", action='store_true', default='false');
#parser.add_argument("-flag_lib", help="directory with meta files for generation of Rosetta cmd-line flag-files", default=__topdir__+"/flag_library" );
#parser.add_argument("-job", help="run-scripts for which job-type? -- default all", choices=job_choices);
#parser.add_argument("-traceback", help="print full traceback in case of error",  action='store_true', default=False )
#parser.add_argument("-label", help="give a special label to this setup -- reflected in file-name for method-options", default='standard' )
#parser.add_argument("-native", help="supply a native pdb for RMSD calculation" )
#keep last!
#parser.add_argument("-method", help="choose algorithm; use together with -h to see additional method-specific options", choices=method_choices );
#parser.add_argument("-database", help="rosetta database location", default='$HOME/rosetta/rosetta_database');
#parser.add_argument("-binaries_prefix", help="rosetta database location", default='$HOME/rosetta/rosetta_source/bin');
#parser.add_argument("-platform", help="which platform", choices=['linux','macos'], default='linux' );
#parser.add_argument("-comp", help="which compiler", choices=['gcc','icc'], default='gcc' );
#parser.add_argument("-extras", help="which extras", choices=['default','static','mpi'], default='default' );
#parser.add_argument("-overwrite", help="overwrite run-directory", default=False, action='store_true' )

# this cqn be run with a pre-fabricated "run_dir" then it is similar to homo_bench setup
# presetup only:  loose files --> target_dir
# runsetup : target_dir + target_options --> run_dir
# full-setup: loose_files --> target_dir + target_options --> run_dir
# target options allows to have target specific stuff: RDC present, exclude residues x1 x2 from scoring, etc.
# for re-setup I can specify a different target_option that is stored in the target_dir (e.g., autoNOE vs ILV-fix )

