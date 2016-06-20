#!/usr/bin/python

import sys
import os

def stripped_file(filename):
        file=open(filename,"r")
        try:
                return ''.join([line.strip() for line in file])
        finally:
                file.close()

original=stripped_file(sys.argv[1])

aa_list = [ '{"'+aa.upper()+'",'+str(i+1)+'}' for i, aa in enumerate(original) ]

parent_table='parent_table={' + ','.join(aa_list) + "}"

print parent_table
