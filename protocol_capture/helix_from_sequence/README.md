Helix from sequence Protocol Capture
====================================
KEYWORDS: MEMBRANES STRUCTURE PREDICTION

Follow the steps below. Example inputs and outputs are provided in the respective folders. 

1) Adjust the paths for the Rosetta executable and database in the cmd_model_tmhelix.sh script

2) Run the script:

   > ./cmd_model_tmhelix.sh

3) Examine the output and compare it to the provided one. It is always recommended to 
   visualize the protein with respect to the membrane bilayer. Example scripts for 
   this are provided in the 006_mp_visualize folder.
   
   Note that in PyMOL line or stick representation, the membrane residue can show up 
   as straight connecting lines between the individual virtual atoms. 
   This artifact can be safely ignored. 
