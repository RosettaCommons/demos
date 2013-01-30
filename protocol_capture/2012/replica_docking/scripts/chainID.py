#!/usr/bin/env python

import subprocess
from os import dup2,path
from os.path import exists
from operator import add
from math import sqrt
from os.path import *
import argparse
import sys
import copy
import shutil
### toolbox library
#import library
from os import mkdir , makedirs

from warnings import *

import traceback

from Bio.PDB.PDBParser import *
from Bio.PDB import PDBIO
from numpy import square

# parser = argparse.ArgumentParser(prog=basename(__file__),
#                                  fromfile_prefix_chars='@',
#                                  description="setup run and target directories",
#                                  add_help=True)
# parser.add_argument("-pdbin", help="pdb file");
# parser.add_argument("-pdbout",help="pdb out")
# parser.add_argument("-id", help="new chain ID",default='B')
# args = parser.parse_args()

parser=PDBParser()
s=parser.get_structure( "t000_", sys.argv[1] )
for chain in s[0]:
    chain.id=sys.argv[2]
io=PDBIO()
io.set_structure( s )
io.save( "%s%c.pdb"%(sys.argv[1],sys.argv[2]) )


