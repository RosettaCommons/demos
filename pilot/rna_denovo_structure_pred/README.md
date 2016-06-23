KEYWORDS: STRUCTURE_PREDICTION NUCLEIC_ACIDS RNA
The goal of this tutorial is to do de novo structure prediction for an
small RNA. The target structure is the SARCIN/RICIN LOOP FROM RAT 28S
R-RNA (PDB 430D). The PDB and it's fasta file are included in the
starting_files directory. 

A good starting point would be the manual:
http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/app_RNA_denovo.html
There are also prep scripts in mini/demo/rna/ which help make the
structure Rosetta friendly. They should probably be included again
here. Other good references are the integration and scientific tests:
mini/test/scientific/tests/rna_denovo
mini/test/integration/tests/rna_denovo

[TROUBLE SHOOTING]

##################################################################
Databases 
##################################################################
minirosetta_database v. rosetta_database 

Someone should double check the public release for 3.3. It appears as though the /rosetta_database/chemical directory is not present 

[!!PUBLIC RELEASE FAILS!!]

##################################################################
Residue Names
##################################################################

The PDB uses capital letters for the residue names in a .fasta file. rna_denovo requires lower case letters



##################################################################
Directory Structure of this Demo
##################################################################


rosetta_inputs/
	run_flags
			-native ./rosetta_inputs/native.pdb
			-fasta ./rosetta_inputs/4d30_.fasta
			-params_file run.prm
			-nstruct 1 [Number of predictions made]
			-out::file::silent 4ds0.out
			-cycles 1000  [30,000 cycles is necessary for larger structures]
			-minimize_rna
			-filter_lores_base_pairs
			-output_lores_silent_file
			-dump [This outputs the PDB of the final predictions]
			-mute core.io.database

	The run_flags file contains command line options

	430D_.fasta
	This is the fasta file that conatins the sequence of the RNA molecule. MAKE SURE YOU USE LOWERCASE RESIDUE NAMES

	native.pdb 
	If hte structure is known it can be used to calculate RMSDs from decoys at the end

	run_params.prm 
	Defines an chain cuts and explicit base pairs

scripts/
	run.sh 
		command used to run rna_denovo


##################################################################
Output created
##################################################################

The output generated in this demo is the 4ds0.out file contianing the score for each prediction genertated and the corresponding PDB files S_XXXXXX.pdb that contain the all atom predictions.



