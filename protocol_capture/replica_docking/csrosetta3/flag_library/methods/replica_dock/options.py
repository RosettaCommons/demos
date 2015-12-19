##-*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t;python-indent:2 -*-'
## make mammoth structure alignments

import string
from glob import glob
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
sub_method_code = flag_lib+"/methods/_docking_base/replica_docking_options.py"
if exists( sub_method_code ):
    exec open(sub_method_code, 'r' )
else:
    print "CANNOT FIND METHOD CODE %s"%sub_method_code
    exit()

# definition of options for the method RASREC



class ReplicaDockingMethod(ReplicaDockingBaseMethod):

	def __init__(self,name,path):
		print "ReplicaDockingMethod: ", name
		ReplicaDockingBaseMethod.__init__(self,name,path)

	def make_target_flags(self, run, setup, filename, flags, subs  ):
		ReplicaDockingBaseMethod.make_target_flags( self, run, setup, filename, flags, subs )

method = ReplicaDockingMethod(method_name, method_path)
