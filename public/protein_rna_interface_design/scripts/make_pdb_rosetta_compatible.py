#!/usr/bin/python

from os import system,popen
from os.path import exists,dirname,basename,expanduser,abspath
from sys import exit, argv
import string
from time import sleep
import glob
import time
import fileinput
import os

input_pdb_filename=argv[1]

assert(exists(input_pdb_filename))

system("cp %s %s "  %(input_pdb_filename, "BACKUP_" + basename(input_pdb_filename) ) )

input_pdb_lines=open( input_pdb_filename ).readlines()

###############################################################
temp_pdb_filename="rna_temp.pdb"

temp_pdb = open( "temp_pdb_filename", 'w' ) #Where the dag commands are ouputted.

for line in input_pdb_lines:
	line_list=line.split()
	if(len(line_list)<=1): continue
	res_name=line_list[3]
	if(res_name=="C" or res_name=="U" or res_name=="A" or res_name=="G"):
		temp_pdb.write(line)

temp_pdb.close()

rna_ready_pdb_filename="rna_temp_RNA_A.pdb"

system("make_rna_rosetta_ready.py %s " %(temp_pdb_filename) ) 

output_pdb_filename=input_pdb_filename.replace(".pdb", "_rosetta_ready.pdb")

output_pdb = open( output_pdb_filename, 'w' ) 

rna_pdb_lines=open( rna_ready_pdb_filename).readlines()

for line in input_pdb_lines:
	line_list=line.split()
	if(len(line_list)<=1): continue
	res_name=line_list[3]
	if((res_name=="C" or res_name=="U" or res_name=="A" or res_name=="G")==False):
		output_pdb.write(line)



for line in rna_pdb_lines:
	line_list=line.split()
	if(len(line_list)<=1): continue
	res_name=line_list[3]
	output_pdb.write(line)

output_pdb.close()


system("renumber_pdb_in_place.py %s " %(output_pdb_filename))


new_rna_pdb_lines=open( output_pdb_filename ).readlines()

output_pdb = open( output_pdb_filename, 'w' ) 

found_TER=False
for n in range(len(new_rna_pdb_lines)):
	line=new_rna_pdb_lines[n]
	line_list=line.split()

	res_name=line_list[3]

	if(found_TER==False):
		if((res_name=="rC" or res_name=="rU" or res_name=="rA" or res_name=="rG")):
			print "n=%d, res_name= %s, line=%s" %(n, res_name, line)
			output_pdb.write("TER\n")
			found_TER=True

	output_pdb.write(line)


output_pdb.close()

