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
from glob import glob
#from sys import argv,stderr,exit
#import sys
from os import popen,system,fdopen,mkdir,makedirs
from os import dup2,path
from os.path import exists
from operator import add
from math import sqrt
from os.path import basename
import argparse
import sys
import shutil
### toolbox library




import traceback

#print "ABRELAX: ", method_path, flag_lib, method_name
sub_method_code = flag_lib+"/methods/_denovo_base/denovo_options.py"
if exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()


# definition of options for the method RASREC
# no extra options necessary
# group=parser.add_argument_group(title='rasrec', description='options for rasrec method')
if 'group' in locals():
    group.add_argument('-cs',help='add chemical shift score to final output',metavar='<file>.tab')
    group.add_argument("-fix_topol", help="provide topology file (obtained with r_pdb2top from pdb-file" );

if 'run_group' in locals():
    run_group.add_argument('-normalize', help='normalize scores by collecting statistics during first rounds of simulation', action='store_true', default=False )
    run_group.add_argument('-quick_relax',help='conduct quick relax during centroid',action='store_true',default=False )
    run_group.add_argument('-pool_size', help='change number of structures in archive -- less means faster calculations', type=int, default=100 )
    run_group.add_argument('-chem_shift_weight', help='set the weight for chemical shift score in filtering [default 5]', default=5, type=float )

class RasrecBaseMethod(DenovoBaseMethod):
    def __init__(self,name,path):
        DenovoBaseMethod.__init__(self,name,path)
        self.option2dir['cs']=self.nmr_sub_dir
        self.option2dir['fix_topol']='structural_knowledge'
        self.SUBS=[]

    def make_target_flags(self, run, setup, filename, flags, subs ):
#        from library import *
        args=self.get_args();
        run.add_subst("CM_AUTO_NSTRUCT", '"-out:nstruct $NSTRUCT "')
        if args.cs:
            flags.write("@flags_cs_rescore\n");
            run.add_subst("CM_CHEM_SHIFTS",setup.cm_path( args.cs ) );

        if args.fix_topol:
            flags.write("#start in stage3 and use this topology only\n")
            flags.write("-iterative:initial_beta_topology %s\n"%
                        setup.cm_path(args.fix_topol))
            flags.write("-iterative:recompute_beta_Naccept 2000\n")
            flags.write("\n");
        DenovoBaseMethod.make_target_flags(self, run, setup, filename, flags, subs)

    def setup_file_library( self ):
        DenovoBaseMethod.setup_file_library( self )
        args=self.get_args()

        fl = self.file_library
        path = flag_lib+"/methods/_rasrec_base/"

        fl.provide_file( "patches", path, "nmr_pool_patch" )
        fl.provide_file( "patches", path, "super_quick_relax.patch" )
#        fl.provide_file( "patches", path, "nmr_relax_patch" )

        fl.provide_file( "flags", path, "flags_iterative" )
        fl.provide_file( "flags", path, "flags_cs_rescore" )
        fl.provide_file( "others", path, "noe_super_quick_relax.txt" )
        fl.add_string( "commandline", "@flags_iterative -run:archive");
        fl.add( "flags", "flags_nmr_patches", "-iterative:fa_score_patch  @@nmr_pool_patch" )
        fl.add( "flags", "flags_nmr_patches", "-iterative:cen_score_patch @@nmr_pool_patch" )

        if args.fix_topol:
            fl.override("flags", "flags_iterative", "-iterative:max_nstruct -1 -1 0 0 0 0 0 0" )

        if args.quick_relax:
            fl.add("flags","flags_iterative", "-iterative:centroid_quickrelax" );

        if args.pool_size:
            fl.override("flags", "flags_iterative", "-iterative:pool_size %d"%args.pool_size )

        if args.normalize:
            fl.add("flags","flags_iterative", "\n#normalize scores:\n-iterative:normalize:activate" );
            fl.add("flags","flags_iterative", "#collect structures for normalization in extra archive\n-iterative:normalize:extra_archive" );
            fl.add("flags","flags_iterative", "#increase this to keep adding structures to the variance-archive\n-iterative:normalize:keep_adding 0.0");
            fl.add("flags","flags_iterative", "-iterative:normalize:nstruct 1000" );
            fl.add("flags","flags_iterative", "#start when all nstruct decoys are collected\n-iterative:normalize:start 0" );
            fl.add("flags","flags_iterative", "#use also for sampling\n-iterative:normalize:sampling" ); #use also for sampling
            fl.add("flags","flags_iterative", "-iterative:normalize:force_zero none" );
            fl.add("flags","flags_iterative", "-iterative:normalize:lower_quartile false" );

        if args.chem_shift_weight:
            fl.override("flags","flags_cs_rescore","-iterative:cenpool_chemicalshift_weight %f"%args.chem_shift_weight)
            fl.override("flags","flags_cs_rescore","-iterative:fapool_chemicalshift_weight %f"%args.chem_shift_weight)
#method = Method(method_name, method_path)
