#!/usr/bin/python

from sys import stdout,argv
from os import system

def stripchain(actualpdbname, out, startchain, newchain, chainspecified):
    lines = open(actualpdbname,'r').readlines()
#    out = open(actualpdbname_newchain,'w')
    for i in range( len(lines)):
        line = lines[i]
        if line.count('ATOM') and (line[21:22] == startchain or not chainspecified):
            line = line[0:21]+newchain+line[22:]
            out.write(line)
    out.close()


actualpdbnames = argv[1:-1]
#actualpdbname_newchain = stdout

for actualpdbname in actualpdbnames:
    out = open('tmp','w')
    newchain = argv[-1]
    chainspecified = 0
    startchain = '_'
#    if len(argv)>3:
#        startchain = argv[-2]
#        if len(startchain) == 1:
#            chainspecified = 1

    if startchain == '-' or startchain =='_':
        startchain = ' '

    if newchain == '-' or newchain =='_':
        newchain = ' '

    stripchain(actualpdbname, out, startchain, newchain, chainspecified)

    command = 'mv tmp '+actualpdbname
    system(command)
