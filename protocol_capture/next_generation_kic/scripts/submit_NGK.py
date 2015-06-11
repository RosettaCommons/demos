#!/usr/bin/python
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1-500
#$ -l arch=lx24-amd64
#$ -l mem_free=5.5G
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
	
home_dir = os.path.expanduser("~")

id = int(sge_task_id)

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

	
bin_path = rosetta_dir + "loopmodel." + bin_ext



if len(sys.argv) < 3:
	print "Error -- you need to provide the list of input PDB files and the target directory core name [will be combined with the PDB name]. The loop file will be derived from the PDB id. You can specify additional flags for KIC after these options (e.g. '-score:weights soft_rep_design')"
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
loop_file = pdb_core + ".loop"
if (not(os.path.isfile(loop_file))):
	print "Error -- loop file %s not found!" % loop_file
	exit(-1)

fa_params = pdb_core + ".fa.params"
cen_params = pdb_core + ".cen.params"

res_dir_core = sys.argv[2]
res_dir = res_dir_core + "_" + pdb_core
if not(os.path.isdir(res_dir)):
	os.popen('mkdir %s' % res_dir) ## make sure the directory exists or create it, otherwise the KIC runs would be for naught
#-

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

outfile_name = pdb_id.split('/')[-1].split(".")[0]+"_committed_"+str(id)+"_"+res_dir_core

local_os = sys.platform
rosetta_cmd_sep = "::"

if (local_os == "darwin"):
	rosetta_cmd_sep = ":"
else:
	home_dir = os.path.expanduser("~") # on the cluster, automatically fetch home directory
#-

args = [
	bin_path,
	"-database %s" % rosetta_db,
	#	"-s %s" % pdb_id, ## KIC doesn't use -s
	"-in:file:fullatom",
	"-loops:loop_file %s" % loop_file,
	"-loops:remodel perturb_kic",
	"-loops:refine refine_kic",
	"-in:file:native %s" % pdb_id,
	"-in:file:s %s " % pdb_id,
	"-out:prefix %s" % outfile_name,
	"-overwrite",
	"-out:pdb_gz",
	"-nstruct 1",
	"-out:path %s" % res_dir,
	"-mute core.io.pdb.file_data",
	"-ex1 -ex2 -extrachi_cutoff 0",
	"-loops:outer_cycles 5",
	"-loops:kic_bump_overlap_factor 0.36",
	"-legacy_kic false",
	"-kic_min_after_repack true",
	"-corrections:score:use_bicubic_interpolation false",
	"-loops:kic_rama2b", 
	"-loops:kic_omega_sampling",
	"-allow_omega_move",
	"-loops:ramp_fa_rep",
	"-loops:ramp_rama",
]

args.extend(addtl_cmds)

print bin_path, args

outfile = open("%s.log"%(outfile_name), "w")
process = subprocess.Popen(args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True)
returncode = process.wait()
outfile.close()

## gzip outfile and move into result directory
gzip = "/usr/bin/gzip"
os.popen(gzip + " " + outfile_name + ".log")
os.popen("mv %s.log.gz %s/" % (outfile_name, res_dir))

## -out:path doesn't seem to be respected at the moment -- move files by hand for now, but fix as soon as possible
os.popen("mv %s%s_0001.pdb.gz %s/" % (outfile_name, pdb_id[0:len(pdb_id)-4], res_dir))

time_end = time.time()

print "Seconds:", time_end - time_start
print "Time:", datetime.timedelta(seconds = time_end - time_start)
print "Summary:", socket.gethostname(), time_end - time_start, datetime.timedelta(seconds = time_end - time_start)
