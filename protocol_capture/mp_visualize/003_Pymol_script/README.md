Pymol script visualization Protocol Capture
===========================================
KEYWORDS: MEMBRANES, ANALYSIS

Visualizing the membrane protein with a PyMOL script. Simply run PyMOL and the script with

   > pymol input/4P79_opm_A.pdb visualize_membrane.pml

The script does not require a span file and (1) either uses the default membrane at the origin
(0, 0, 0) with the normal in the z-dimension or (2) uses the information stored in the membrane 
residue MEM. The script requires 'membrane_planes.py' in the same directory. 