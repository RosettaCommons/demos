#!/usr/bin/env python2.7

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
import traceback

import argparse
import sys

# Zhe, Oct 6,2011, histogram

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog=basename(__file__),
                                 fromfile_prefix_chars='@',
                                 description='histogram of given file of given parameter',
                                 add_help=True)

parser.add_argument("-files",nargs='*', help = ".out files")
parser.add_argument("-var", help = "histogram of which variable",default = 'rms')
parser.add_argument("-bin", help = "bins of the histogram", default = '20')
#parser.add_argument("-colors",nargs='*',  help = "give a color for each file", default = 'blue')
parser.add_argument("-title",help = "title of the figure")
parser.add_argument("-loc",help = "location of legend",default = 0)

args=parser.parse_args()

def readFile(file,var):
    var_idx = -1
    var_data = []
    lines = open(file,'r').readlines()
    for line in lines:
        tags = line.split()
        if var_idx < 0:
            if tags[0]=='SCORE:' and var in tags:
                var_idx = tags.index(var)
        else:
            if tags[0]=='SCORE:' and tags[var_idx]!=args.var:
                var_data.append(float(tags[var_idx]))
    return var_data

data = []
for i in range(len(args.files)):
    data.append( readFile(args.files[i],args.var) )
    
#n,bins,patches = plt.hist(data,int(args.bin),normed=1,color=args.colors,label=args.files)
n,bins,patches = plt.hist(data,int(args.bin),normed=1,label=args.files)
# print n
# print bins

# plt.plot(bins[:-1]+(bins[1]-bins[0])/2.0,n,'r--')
plt.legend(loc=args.loc)
#plt.grid(True)
plt.xlabel(args.var)
plt.title(args.title)
plt.show()




