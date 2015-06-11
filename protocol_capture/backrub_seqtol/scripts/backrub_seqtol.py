#!/usr/bin/python
# These are Sun Grid Engine options used on the UCSF cluster, some may need to
# be changed for your local environment
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=10:00:00
#$ -l arch=lx24-amd64

import socket
import sys

# Print the python/host information as soon as possible
print "Python:", sys.version
print "Host:", socket.gethostname()

import datetime
import os
import subprocess
import time

# Users should change these variables to match their local environment
home_dir = os.path.expanduser("~")
no_hb_env_dep_weights_file = os.path.join(os.getcwd(), "input_files", "standard_NO_HB_ENV_DEP.wts")
backrub_bin = os.path.join(home_dir, "rosetta/rosetta-3.2/rosetta_source/bin/backrub.linuxgccrelease")
seqtol_bin = os.path.join(home_dir, "rosetta/rosetta-3.2/rosetta_source/bin/sequence_tolerance.linuxgccrelease")
database_dir = os.path.join(home_dir, "rosetta/rosetta-3.2/rosetta_database")

# Score function modifications (Smith & Kortemme 2010 set all these to True)
no_hb_env_dep = False
no_score12_patch = False
histidine_offset = True

# Settings dependent on score function modifications
score_weights = "standard"
if no_hb_env_dep:
	score_weights = no_hb_env_dep_weights_file
score_patch = "score12"
if no_score12_patch:
	score_patch = ""
ref_offsets = []
if histidine_offset:
	ref_offsets = ["HIS", "1.2"]

# Determine the starting PDB structure from either the PDB_PATH environment
# variable, or the first command line argument
if os.environ.has_key("PDB_PATH"):
	pdb_path = os.path.abspath(os.environ["PDB_PATH"])
elif len(sys.argv) >= 2:
	pdb_path = os.path.abspath(sys.argv[1])
else:
	sys.exit("No start PDB file specified")

# Determine the iteration from either the SGE_TASK_ID environment variable or 
# the second command line argument
iteration = 1
if os.environ.has_key("SGE_TASK_ID"):
	iteration = int(os.environ["SGE_TASK_ID"])
elif len(sys.argv) >= 3:
	iteration = int(sys.argv[2])

# Check if the start PDB has valid printf syntax and select the correct PDB
index = 1
if pdb_path.find("%") >= 0 and os.path.exists(pdb_path % index):
	pdb_paths = [pdb_path % index]
	index += 1
	while os.path.exists(pdb_path % index):
		pdb_paths.append(pdb_path % index)
		index += 1
	pdb_path = pdb_paths[(iteration-1) % len(pdb_paths)]

# Determine the PDB chain identifiers
pdb_file = open(pdb_path)
pdb_lines = pdb_file.readlines()
pdb_file.close()
atom_lines = [line for line in pdb_lines if line.find("ATOM  ") == 0 and len(line) >= 22]
pdb_chains = []
for line in atom_lines:
	if not line[21] in pdb_chains:
		pdb_chains.append(line[21])

# Set the fitness function master weights based on the number of chains
num_chains = len(pdb_chains)
fitness_function_weights = ["1"] + ["1"]*num_chains + ["2"]*(num_chains*(num_chains-1)/2)
# Customize your weights by uncommenting this line
#fitness_function_weights = ["your", "weights", "go", "here"]

# Get the base path from either the PDB file, the INPUT_PATH environment
# variable, or the third command line argument
input_path = os.path.splitext(pdb_path)[0]
if os.environ.has_key("INPUT_PATH"):
	input_path = os.path.abspath(os.environ["INPUT_PATH"])
elif len(sys.argv) >= 4:
	input_path = os.path.abspath(sys.argv[3])

# Determine input/output directories and filenames from the starting PDB
pdb_name = os.path.split(os.path.splitext(pdb_path)[0])[-1]
input_name = os.path.split(input_path)[-1]
output_dir = os.path.join(os.getcwd(), "output", input_name)
output_prefix = input_name + "_%04i" % (iteration)

# Print job summary information
print "Iteration:", iteration
print "Database:", database_dir
print "PDB:", pdb_path
print "PDB Chains:", ", ".join(['"'+chain+'"' for chain in pdb_chains])
print "Fitness Function Weights:", ", ".join(fitness_function_weights)
print "Inputs:", input_path
print "Output Directory:", output_dir
print "Output Prefix:", output_prefix
print "No hydrogen bond environment dependence:", no_hb_env_dep
print "No score12 Patch:", no_score12_patch
print "Histidine Offset:", histidine_offset

backrub_args = [
	backrub_bin,
	"-database", database_dir,
	"-s", pdb_path,
	"-ignore_unrecognized_res",
	"-ex1", "-ex2", "-extrachi_cutoff", "0",
	"-out:prefix", output_prefix,
	"-mute", "core.io.pdb.file_data",
	"-backrub:ntrials", "10000",
	"-score:weights", score_weights,
	"-score:patch", score_patch,
]

# Add backrub resfile to the command line if it exists
backrub_resfile = os.path.join(input_path + "_backrub.resfile")
if os.path.exists(backrub_resfile):
	backrub_args += ["-resfile", backrub_resfile]
	print "Backrub Resfile:", backrub_resfile

# Add minimize movemap to the command line if it exists
minimize_movemap = os.path.join(input_path + "_minimize.movemap")
if os.path.exists(minimize_movemap):
	backrub_args += ["-backrub:minimize_movemap", minimize_movemap]
	print "Minimize Movemap:", minimize_movemap

# Add perturb movemap to the command line if it exists
perburb_movemap = os.path.join(input_path + "_perturb.movemap")
if os.path.exists(perburb_movemap):
	backrub_args += ["-in:file:movemap", perburb_movemap, "-sm_prob", "0.1"]
	print "Perturb Movemap:", perburb_movemap

seqtol_args = [
	seqtol_bin,
	"-database", database_dir,
	"-s", output_prefix + "_low.pdb",
	"-ex1", "-ex2", "-extrachi_cutoff", "0",
	"-seq_tol:fitness_master_weights"] + fitness_function_weights + [
	"-ms:generations", "5",
	"-ms:pop_size", "2000",
	"-ms:pop_from_ss", "1",
	"-ms:checkpoint:prefix", output_prefix,
	"-ms:checkpoint:interval", "200",
	"-ms:checkpoint:gz",
	"-ms:numresults", "0",
	"-out:prefix", output_prefix,
	"-score:weights", score_weights,
	"-score:patch", score_patch,
]

# Add ref_offsets to commane line if they are defined
if len(ref_offsets):
	seqtol_args += ["-score:ref_offsets"] + ref_offsets

# Add resfile to command line if it exists
seqtol_resfile = os.path.join(input_path + "_seqtol.resfile")
if os.path.exists(seqtol_resfile):
	seqtol_args += ["-resfile", seqtol_resfile]
	print "Sequence Tolerance Resfile:", seqtol_resfile

# Make the output directory and make it the working directory
try:
	os.makedirs(output_dir)
except:
	if not os.path.exists(output_dir):
		print "Couldn't make output directory!"
		sys.exit(1)
os.chdir(output_dir)

time_start = time.time()

# Run the backrub app if it hasn't already been run
if not (os.path.exists(output_prefix + "_low.pdb") or os.path.exists(output_prefix + "_low.pdb.gz")):
	print " ".join(backrub_args)
	sys.stdout.flush()

	outfile = open(output_prefix + "_backrub.out", "w", 0)
	process = subprocess.Popen(backrub_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds=True)
	returncode = process.wait()
	outfile.close()

	# Remove the "last" structure and rename the "low" structure"
	os.remove(output_prefix + pdb_name + "_0001_last.pdb")
	os.rename(output_prefix + pdb_name + "_0001_low.pdb", output_prefix + "_low.pdb")
else:
	print "Skipping backrub (Already Completed)"

# Print backrub running time statistics
time_backrub_end = time.time()
print "Time (backrub):", datetime.timedelta(seconds = time_backrub_end - time_start)

# Run the sequence_tolerance app if the PDB file hasn't been compressed
if os.path.exists(output_prefix + "_low.pdb"):
	print " ".join(seqtol_args)
	sys.stdout.flush()

	seqtol_output_filename = output_prefix + "_seqtol.out"
	output_index = 0
	# Find an output filename that hasn't been used
	while os.path.exists(seqtol_output_filename):
		output_index += 1
		seqtol_output_filename = output_prefix + "_seqtol%i.out" % output_index
	outfile = open(seqtol_output_filename, "w", 0)
	process = subprocess.Popen(seqtol_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds=True)
	returncode = process.wait()
	outfile.close()

	# Compress the output structure
	process = subprocess.Popen(["/bin/gzip", output_prefix + "_low.pdb"], close_fds=True)
	returncode = process.wait()
else:
	print "Skipping sequence_tolerance (Already Completed)"

# Print running time statistics
time_end = time.time()
print "Time (seqtol):", datetime.timedelta(seconds = time_end - time_backrub_end)
print "Seconds:", time_end - time_start
print "Time:", datetime.timedelta(seconds = time_end - time_start)
print "Summary:", socket.gethostname(), time_end - time_start, datetime.timedelta(seconds = time_end - time_start)
