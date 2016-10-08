MP Transform Optimize Protocol Capture
======================================
KEYWORDS: MEMBRANES ANALYSIS


Follow the steps below. Example inputs and outputs are provided in the respective folders. 
Three different options are provided for variable flexibility:

   cmd_mptransform.sh           transforms the protein into the membrane
   cmd_mptransform_optimize.sh  transforms the protein into the membrane and optimizes the 
	                             embedding with the high-resolution membrane score function

1) Adjust the paths for the Rosetta executable and database in the cmd_mptransform<...>.sh scripts.
   Make sure you understand the options in the cmd files. 

2) Run the scripts:

   > cmd_mptransform.sh
   > cmd_mptransform_optimize.sh

3) Check your output: visualize the positions of the input and output structures in PyMOL with the
   visualize_membrane.pml script using the command line
   
   > pymol <output>.pdb visualize_membrane.pml 

