#!/usr/bin/python
## make mammoth structure alignments

import string
from glob import glob
from sys import argv,stderr,exit,stdout
from os import popen,system
from os.path import exists
from operator import add
from math import sqrt

#############################
def Help():
    print '\n'
    print '-'*75
    print 'USAGE: %s fa_cst'%argv[0]
    print '\n will map all QX to CB and add 1.5 padding per methyl .. result printed to screen'
    print '\n\n\n\n'
    exit()

#################
def rename_atom( atom ):
    if atom=="HA2": 
	return "CA",2
    if atom=="HA3":
	return "CA",2
    if atom[0:2]=='HA':
	return 'CA',2	
    if atom[0:2]=='QQ':
	return 'CB',6
    if atom[0:2]=='QA':
	return 'CA',2	
    if atom[0]=='Q':
        return 'CB',3 
    if atom[0:3]=='HE1':
        return 'CB',1
    if atom[0:2]=='HG':
        return 'CB',1
    if atom[0:2]=='HB':
        return 'CB',1
    if atom[0:2]=='HH':
	return 'CB',1
    if atom[0:2]=='HE':
	return 'CB',1
    if atom[0:2]=='HD':
	return 'CB',1
    if atom[0:2]=='HZ':
	return 'CB',1
    if atom[1:3]=='HZ':
	return 'CB',1
    return atom, 1

if len(argv) <=1:
    Help()




file = argv[1]

lines = open(file).readlines()

for line in lines:
    if (line[0]=="#"): continue;
    l=string.split(line);
    try:
        resi1 = int( l[2] );
        resi2 = int( l[4] );
    except:
        print '\n\n expected integer in col 3 and 5 of:\n'
        print l[2],l[4],'\n',l
        exit()

    atom1,met1 = rename_atom( l[1] );
    atom2,met2 = rename_atom( l[3] );
    
    bound_low = float( l[6] );
    try:
        bound_up = float( l[7] );
        bound_sd = float( l[8] );
    except:
        print '\n\n expected float in col 8 and 9 of:\n'
        print l[7],'\n',l
        exit()

    pad = met1 * met2
    hyd = (met1 > 1) + (met2 > 1 )	
    #print line, pad, hyd, pad**(1.0/6)
    bound_sd_ori = bound_sd
    bound_up_ori = bound_up
    bound_sd = bound_sd + 0.5*((pad-1)**(1.0/6))
    bound_up = (bound_up+0.5*hyd)*(pad**(1.0/6))

    print 'AtomPair %5s %5d %5s %5d BOUNDED %5.3f %5.3f %5.3f NOE; rawdata %5.3f %5.3f'%(atom1,resi1,atom2,resi2,bound_low,bound_up,bound_sd,bound_up_ori, bound_sd_ori)







