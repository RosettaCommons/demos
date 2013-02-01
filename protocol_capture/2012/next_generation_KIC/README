### General Information ##################
# Your name:
Amelie Stein, Tanja Kortemme

# Protocol Name:
Next-generation KIC

# Publication Describing the Method:
Stein A, Kortemme T (2013)
Improvements to robotics-inspired conformational sampling in Rosetta.
PLoS ONE (submitted)


# Brief Description:

### setting up the demo: #######

Before running the demo, make sure to change the following variables in your local environment:
PATH_TO_EXE     ## path to directory with Rosetta binaries
ROSETTA_BINARY  ## extension of Rosetta binaries, e.g. linuxgccrelease
PATH_TO_DB 	## path to Rosetta database

### running #########
# Important Flags, other flags that can be used, required flag, brief description of each flag or option:
# These are all found in input_files/flags
#io flags:
-s 1a8d_MinPacked.pdb                     	      	# The starting structure -- must have residues for the segment to be remodeled, but these don't need to have meaningful coordinates
-loops:loop_file 1a8d.loop 	  	      		# definition of the loop to be remodeled, file format description at http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d1/d49/loopmodeling.html
-out:pdb_gz                               	      	# compress output structures
-in:file:native 1a8d_MinPacked.pdb	  	      	# native or reference structure for RMSD calculation -- if this flag is not set, an RMSD of 0 will be reported

#number of structures to produce                        
#for demo:
-nstruct 1                                	      	# number of structures to produce 
#for production runs, generate at least 500 structures, depending on the cluster setup in independent simulations with -nstruct 1 each or with -nstruct 500

#general kinematic closure loop modeling flags:
-in:file:fullatom
-loops:remodel perturb_kic
-loops:refine refine_kic
#-run:test_cycles  			  	      	# fast execution, uncomment for testing purposes only
-loops:outer_cycles 5			  	      	# for production runs
-kic_bump_overlap_factor 0.36		  	      	# reduces the threshold for permitted clashes in initial loop closures (default value is 0.4)
-legacy_kic false			  	      	# remove a slight bias in the initial KIC implementation towards sampling the C-termial part of the remodeled segment more frequently
-kic_min_after_repack true		  	   	# minimize rotamers after repacking (happens every -loops:repack_period iterations after loop closure, default 20)
-corrections:score:use_bicubic_interpolation false	# do not use bicubic interpolation, it has been observed to adversely affect the fraction of sub-Angstrom conformations on the 12-residue benchmark set

#next-generation KIC flags as described in Stein & Kortemme, PLoS ONE, 2013
-loops:kic_rama2b					# note: this requires 5G of memory
-loops:kic_omega_sampling
-allow_omega_move true
-loops:ramp_fa_rep
-loops:ramp_rama

#packing flags
-ex1
-ex2
-extrachi_cutoff 0


# Example Rosetta Command Line:
$PATH_TO_EXE/loopmodel.$ROSETTA_BINARY -database $PATH_TO_DB @flags

# Overall protocol execution (demo)

1. scripts/pre_min_pack.py <list_of_native_structures> <output_keyword> (prepacking step)
This will create a pre-min-packed structure to start the simulation from, as well as a log of the run.
The input is a list of starting structures to be processed, as well as a keyword for the output structure
and log names.
Note that, while standard repacking usually takes a few minutes, with the -min_pack option it can take
an hour or more, depending on the size of the structures. This option was developed by Andrew Leaver-Fay
and will be discussed in a separate publication (Leaver-Fay et al., Meth Enzym, 2013).


2. scripts/submit_NGK.py <input_list> <output_keyword> (loop remodeling step)
This will generate structures in which the selected segment is remodeled (1 for the demo,
use at least 500 for real life problems, more for longer segments). The input is a list
of (pre-packed) structures and a keyword for output naming.
The Rosetta syntax for .loop files is explained at 
http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d1/d49/loopmodeling.html

Note that if a native or reference structure is provided, the backbone RMSD of the remodeled segment
(after superimposition of the fixed parts of the structure) will be reported by loopmodel.release and
can thus be parsed from the standard output, if that is redirected into a file. The script will
perform redirection into a .log file which later is gzipped.

If submitted via qsub on an SGE cluster system, the script will distribute the number of
simulations specified with the -t flag across all structures in the list, i.e., to 
generate 500 structures each for a list of 10 input structures, use -t 1-5000.

An NGK trajectory (generating one model) typically take 15-30min for 12-residue loops.

As noted above, using -loops:kic_rama2b requires around 5GB of memory, specified here
with the -mem_free option for the SGE cluster system. On systems with less memory, the Rama2b flag
may be left out; memory requirements for the other NGK flags are 1GB. However, this will lead
to a decrease in the percentage of sub-Angstrom conformations generated across a simulation, as
described in Stein & Kortemme (2013, Fig. 3B). It is thus recommended to generate more structures 
(e.g. 1000) when leaving out the -loops:kic_rama2b flag.



3. scripts/parse_RMSDs.pl <directory with log files> (parsing step)
This script will extract the scores and RMSDs for each final model from the respective .log files,
and optionally generate Rosetta-energy-vs-RMSD plots. Note that relevant RMSDs will only be reported
if a native or reference structure is provided to the NGK run via the -in:file:native flag.
If no reference structure is available, users should consider clustering of their results to
identify the most commonly sampled conformations (see Figure S1 of Stein & Kortemme, PLoS ONE, 2013).


### versions #########
Latest version applies to svn revision 51851 (Dec 2012)

# Version for other codes used, version for interfaces or libs used (if relevant)

# References to published works using this protocol
Stein A, Kortemme T (2013)
Improvements to robotics-inspired conformational sampling in Rosetta.
PLoS ONE (submitted)


# Other Comments: 



