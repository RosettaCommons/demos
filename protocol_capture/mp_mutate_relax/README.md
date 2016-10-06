MP Mutate Relax Protocol Capture
================================
KEYWORDS: MEMBRANES STRUCTURE_PREDICTION STABILITY_IMPROVEMENT


Follow the steps below. Example inputs and outputs are provided in the respective folders. 
Three different options are provided for variable flexibility:

   cmd_mutate.sh         only mutates and repacks the mutated residue
   cmd_mutate_repack.sh  mutates and repacks within a radius of the mutation site
   cmd_mutate_relax.sh   mutates and refines the entire protein with sidechain optimization. 

1) Adjust the paths for the Rosetta executable and database in the cmd_mutate<...>.sh scripts.
   Make sure you understand the options in the cmd files. 
	
   While mutate and mutate_repack only require a single output model, mutate_relax requires more
   models. Choose a number between 50 and 1,000 for production runs, depending on the computational resources. 

2) Run the scripts:

   > ./cmd_mutate.sh
   > ./cmd_mutate_repack.sh
   > ./cmd_mutate_relax.sh

3) Check your output: the relevant output files have the desired mutation in the filename.
   Scores can be found in the PDBs after the ATOM lines. 

