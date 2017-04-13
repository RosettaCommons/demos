MP Score Protocol Capture
=========================
KEYWORDS: MEMBRANES ANALYSIS

Follow the steps below. Example inputs and outputs are provided in the respective folders. 
Three different options are provided:

   cmd_mpscore_as-is.sh               scores the protein as-is with respect to the fixed membrane
   cmd_mpscore_transform.sh           transforms the protein first into the membrane before scoring
   cmd_mpscore_transform_optimize.sh  transforms the protein first into the membrane, optimizes
                                      the embedding with the high-resolution score function and 
                                      then scores the protein

1) Adjust the paths for the Rosetta executable and database in the cmd_mpscore<...>.sh scripts.
   Make sure you understand the options in the cmd files. 

2) Run the scripts:

   > cmd_mpscore_as-is.sh
   > cmd_mpscore_transform.sh
   > cmd_mpscore_transform_optimize.sh

3) Check your output: visualize the positions of the input and output structures in PyMOL with the
   visualize_membrane.pml script:
   
   pymol <input>.pdb <output>.pdb visualize_membrane.pml
   
   The score files are given in the respective output folders. 
	
4) You can visualize clashes via:

   > ./visualize_clashes.py -p output_score_transform_optimize/4P79_A_0001_transform_optimize.pdb -r
	
   The -r option visualizes the repulsive score term. The script writes a PyMOL script
   visualize_clashes.pml that can be run inside PyMOL after loading the structure file. To do
   this, run the following in a command line terminal:
	
   > pymol output_score_transform/4P79_A_0001_transform.pdb visualize_clashes.pml

