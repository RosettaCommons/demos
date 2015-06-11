#!/usr/bin/python

from sys import argv,exit
from os import popen, system
from os.path import basename
import string



def Help():
    print
    print 'Usage: '+argv[0]+' <silent out file 1> < tag file > '

    exit()


if len(argv)<2:
    Help()


infile = argv[1]
tagfile = argv[2]

tags=[]
f = open( tagfile )
for line in f:
    tt = string.split(line[:-1]);
    tags.append( tt[0])  #gets rid of newline

#icommand = 'head -n 2 '+infile
#system(command)

count = 1
fid = open( infile )
line = fid.readline()
writeout = 0
head_score = ''
head_sequence = ''
while line:
        cols = string.split(line)
	if (len(cols)>1 and cols[0]=='SCORE:'):
            try:
                index = tags.index( cols[-1] )
                writeout = 1
                if ( len(head_sequence) > 1 ):
			print head_sequence
		        head_sequence=''
                if ( len(head_score) > 1 ):
			print head_score
			head_score='' 
	    except ValueError:
                writeout = 0
        if cols[1]=='score': 
	    writeout = 0
	    head_score = line[:-1]
        if cols[0]=='SEQUENCE:':
            writeout = 0
	    head_sequence = line[:-1]        
        if writeout:
            print line[:-1]
        line = fid.readline()



