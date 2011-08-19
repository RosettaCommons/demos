#!/usr/bin/python

from sys import argv,stdout
from os import popen,system
from os.path import exists,basename
from get_sequence import get_sequence
removechain = 0
if argv.count('-nochain'):
    removechain = 1

pdbnames = argv[1:]

for pdbname in pdbnames:
    fasta_line = get_sequence( pdbname, removechain )

    fastaid = stdout
    fastaid.write( '>%s %d\n' % (basename(pdbname), len( fasta_line ) ) )
    fastaid.write( fasta_line )
    fastaid.write( '\n' )
