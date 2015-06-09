### General Information ##################
# Your name:Elizabeth Dong

# Protocol Name:Loop Building for Membrane Proteins with Toplogy Broker

# Brief Description:The Topology Broker application in Rosetta can be used to build loops and flexible regions in membrane proteins. The protocol allows for the inclusion of membrane span and lipophilicity information as well as membrane weight set for the scoring function. This protocol rebuilds the intracellular helical region between TM 5 & 6 and C-terminus for Squid Rhodopsin (PDB ID: 2Z73A).


# if this is a protocol attached to a rosettaCon talk provide:

  # Title, Authors & Lab , Year, Session and Day of talk
  Elizabeth Dong, Anette Schreiber, Karen Gregory, Kristian Kauffman, Jeff Conn, Jens Meiler
  Meiler Lab, Vanderbilt University, Nashville, TN
  08/05/10, Session 9
  
  # Abstract
	  Selective modulators of metabotropic glutamate receptor subtype 5 (mGluR5), a class C G-protein coupled receptor, provide novel treatment strategies for disorders that disrupt cognitive function. Identifying the specific residues on mGluR5 that contact these small molecules would allow for a deeper understanding of the binding interaction and aid in the development of therapeutic compounds. Construction of the mGluR5 model entailed identification of TM segments in the sequence using JUFO9D, modeling the 7 TM helices based on the three mammalian GPCR crystal structures using Rosetta and modeling the loops using the Toplogy Broker application in Rosetta, which allows for the inclusion of membrane span and lipophilicity information as well as the use of a membrane weight set to apply to the scoring function. Residues of mGluR5 critical for the binding of allosteric modulators were determined through Rosetta Ligand docking studies informed by experimental functional data. The experimentally validated models demonstrate the success of Rosetta to model GPCRs.
 
### running #########
# Important Flags, other flags that can be used, required flag, brief description of each flag or option:

-run:protocol broker  #initiate call to broker
-broker:setup ./input_files/setup_broker.tpb #defines constraints on folding protocol
-frag3 ./input_files/aa2Z73A03_05.200_v1_3 #fragment files for folding of protein
-frag9 ./input_files/aa2Z73A09_05.200_v1_3 #fragment files for folding of protein
#Patches to the scoring function ensure that membrane potentials are used in the folding protocol
#make sure to have these either in the local directory or your database directory under scoring/weights
-stage2_patch score_membrane_s2.wts_patch
-stage3a_patch score_membrane_s3a.wts_patch
-stage3b_patch score_membrane_s3b.wts_patch
-stage4_patch score_membrane_s4.wts_patch
#allows setup of membrane options
-abinitio
	-membrane
#options for membrane scoring functions
-membrane
	-no_interpolate_Mpair
	-Menv_penalties
#tells folding protocol to close loops
-close_loops
-non_ideal_loop_closing
#not sure what this does but Yeifan included it.
-score
	-find_neighbors_3dgrid
-no_prof_info_in_silentout #no time-columns appears in score/silent - files
#input files these options should actually be supplied in the command line
-in
	-file
		-fasta ./input_files/2Z73A.fasta
		-spanfile ./input_files/2Z73A.span #generate from Octopus prediction (http://octopus.cbr.su.se/) using /TopologyBroker_GPCR/scripts/octopus2span.pl
		-lipofile ./input_files/2Z73A.lips4 #generate using /TopologyBroker_GPCR/scripts/run_lips.pl
#-out
#	-file
		#output options the results.silent_binary_out should be supplied on command line
		#-silent results.silent_binary_out
		#-silent_struct_type binary
#number of structures to generate supply on command line
#-nstruct 1

# Example Rosetta Command Line:
r_broker.linuxgccrelease -database ~/minirosetta_database -out:file:residue_type_set centroid -out:file:silent rbroker_run1.out -nstruct 1 @flags.txt

# Example Overall Command Line (if overall protocol is run via a script or other program)
n/a

### versions #########
# If checked into trunk, svn revision number: 37327

# If NOT checked into trunk, path and svn revision number: n/a

# Version for other codes used, version for interfaces or libs used (if relevant): n/a

# References to published works using this protocol (you may include submitted, and in press works): n/a

# Other Comments: 
Before running the example, put all score_membrane*.wts_patch in the local directory

To generate *.span file: 
generate from Octopus prediction (http://octopus.cbr.su.se/) using /TopologyBroker_GPCR/scripts/octopus2span.pl

To generate *.lips4 file:
run the script /TopologyBroker_GPCR/scripts/run_lips.pl with the following command line:
run_lips.pl <fasta file> <span file> <path to blastpgp> <path to nr database> <path to alignblast.pl script>
example: run_lips.pl 2Z73A.fasta 2Z73A.span /sb/meiler/Linux2/x86/blast/blast-2.2.18/bin/blastpgp /sb/meiler/scripts/sequence_analysis/db/nr /mini/src/apps/public/membrane_abinitio/alignblast.pl 

Topology Broker only runs in centroid mode as of now. To extract pdb from *.out, use:
~/extract_pdbs.linuxgccrelease -database ~/minirosetta_database/ -in:file:silent rrbroker_run1.out -in:file:residue_type_set centroid -out:file:residue_type_set centroid -out:output