#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1-45
#$ -l arch=linux-x64
#$ -l mem_free=2G
#$ -l netapp=1G

# -t specifies the total number of simulations, which will be distributed across the input PDB list

import socket
import sys

print "Python:", sys.version
print "Host:", socket.gethostname()

import datetime
import os
import subprocess
import time
import sys
import commands

time_start = time.time()

sge_task_id = 1
if os.environ.has_key("SGE_TASK_ID"):
	sge_task_id = os.environ["SGE_TASK_ID"]

id = int(sge_task_id)
	
home_dir = os.path.expanduser("~")


## TODO: read in Rosetta-specific variables as outlined in README
rosetta_dir = home_dir + "/rosetta/rosetta_source/bin/"
if os.environ.has_key("PATH_TO_EXE"):     ## path to directory with Rosetta binaries
	rosetta_dir = os.environ["PATH_TO_EXE"]
#-
bin_ext = "linuxgccrelease"
if os.environ.has_key("ROSETTA_BINARY"):      ## extension of Rosetta binaries, e.g. linuxgccrelease
	bin_ext = os.environ["ROSETTA_BINARY"]
#-
rosetta_db = home_dir + "/rosetta/rosetta_database"
if os.environ.has_key("PATH_TO_DB"):     ## path to Rosetta database
	rosetta_db = os.environ["PATH_TO_DB"]
#-

	
bin_path = rosetta_dir + "fixbb." + bin_ext


if len(sys.argv) < 3:
	print "Error -- you need to provide the list of input PDB files and the outfile keyword (e.g. my_min_packed). You can specify additional flags for fixbb after these options"
	exit(-1)
#-

pdb_lst_file = sys.argv[1]
pdb_lst = []
read_pdbs = open(pdb_lst_file)
for l in read_pdbs:
        pdb_lst.append(l.strip('\n'))
#-
read_pdbs.close()
pdb_id = pdb_lst[id % len(pdb_lst)]
pdb_core = pdb_id.split("_")[0] ## 1cb0_min.pdb -> 1cb0

fa_params = pdb_core + ".fa.params"
cen_params = pdb_core + ".cen.params"

outfile_key = sys.argv[2]

addtl_cmds = []
for i in range(3, len(sys.argv)):
	addtl_cmds.append(sys.argv[i])
#-

if os.path.isfile(fa_params):
	addtl_cmds.extend(["-extra_res_fa", fa_params])
#-

if os.path.isfile(cen_params):
	addtl_cmds.extend(["-extra_res_cen", cen_params])
#-

outfile_name = pdb_id.split('/')[-1].split(".")[0]+"_" + outfile_key + "_"+str(id)

local_os = sys.platform
rosetta_cmd_sep = "::"

if (local_os == "darwin"):
	rosetta_cmd_sep = ":"
#-

#lig_path = database_path+"chemical/residue_type_sets/fa_standard/residue_types/metal_ions/MG.params" ## only required for non-built-in types

args = [
	bin_path,
	"-database %s" % rosetta_db,
	"-in:file:fullatom",
	"-in:file:s %s " % pdb_id,
	"-out:suffix %s" % outfile_key,
	"-min_pack",
	"-packing:repack_only", 
	"-ignore_unrecognized_res", ## to ignore ligands, waters etc. that may be in native PDB structures
	"-overwrite",
	#"-out:pdb_gz",
	"-nstruct 1",
	#"-mute core.io.pdb.file_data",
	"-ex1 -ex2 -extrachi_cutoff 0",
]

args.extend(addtl_cmds)

print bin_path, args

outfile = open("%s.log"%(outfile_name), "w")
process = subprocess.Popen(args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True)
returncode = process.wait()
outfile.close()

## gzip outfile and move into result directory
gzip = "/usr/bin/gzip"
os.popen(gzip + " -f " + outfile_name + ".log")

time_end = time.time()

print "Seconds:", time_end - time_start
print "Time:", datetime.timedelta(seconds = time_end - time_start)
print "Summary:", socket.gethostname(), time_end - time_start, datetime.timedelta(seconds = time_end - time_start)
