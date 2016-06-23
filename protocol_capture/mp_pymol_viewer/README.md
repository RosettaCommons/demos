Visualizing Membranes in PyMOL
==============================

KEYWORDS: MEMBRANES UTILITIES

Author: Rebecca F. Alford (rfalford12@gmail.com)  
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: January 2015  
Rosetta Revision #58069 
PyMOL Version: 1.7.4  

---

Visualizing the position and orientation of a biomolecule in the membrane bilayer
(embedding) is an important part of analyzing membrane models from Rosetta. This visualization 
allows the user to distingish conformations with poor embeddings from conformations with 
native-like/reasonable embeddings. 

This application displays a set of two parallel planes in PyMOL, representing the position 
and geometry of the membrane bilayer. During a simulation, Rosetta extracts the membrane center, 
normal and thickness from the MEM residue, calculates the position of the planes, and sends
this information to PyMOL in real time. PyMOL then uses this information to draw two CGO plane objects
representing the membrane. 

This tool is part of a standalone application and can also be used in combination with any JD2-supported
Rosetta application. 

Publication describing the method: 

* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ 
  (2015) An integrated framework advancing membrane protein modeling and design 
  PLoS ONE (in preparation) 

Executable/Script
-----------------

    Rosetta/main/source/bin/view_membrane_protein.linuxgccrelease

"or"

Pass the `-show_simulation_in_pymol 0` flag with any Rosetta Membrane Framework application

Generating Inputs
-----------------

Two inputs are required for using the standalone visualization app: 

1. A PDB to view 
2. Span file describing the location of trans-membrane spans

Steps for generating these inputs are found below. A set of example inputs can 
also be found in example_inputs/. Here, 1c3w is used as an example: 

1. PDB File: If an output model is not already available from Rosetta, 
   generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2014. An example commandline using 
   1c3w is also provided here: 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/db -in:file:s example_inputs/1c3w_tr.pdb

   For this example, this command will produce 1 output files: 
     = 1c3w_tr.span: Spanfile containing predicted trans-membrane spans

Steps of the protocol
---------------------

Here, we describe the steps required to run the MP_PyMOLViewer protocol. As an example, all steps 
use the PDB 1c3w: 

1.  Required Options: Options (flags) needed to run this application. A file with these flags, pymol_flags, 
    is also provided for 1c3w in this demo: 

        flags                                  descriptions
        --------------------------------------------------------------------------------------------------
        -in:file:s <pdbfile>                   Input PDB Structure: PDB file for protein structure
        -membrane_new:setup:spanfiles          Spanfile describing spanning topology of starting structure 
                                               for full symmetric structure
        -show_simulation_in_pymol 0			       Use the PyMOL viewer to visualize membrane planes for structures
        -keep_pymol_simulation_history 1       Keep pymol frames for making movies/replaying simulations (optional)

2.  Startup the PyMOL PyRosetta Session: 

    1. Open a new session of PyMOL
    2. Run the PyMOLPyrosettaServer.py script using the following command line in the pymol window:  

            run /path/to/Rosetta/main/source/src/python/bindings/PyMOLPyRosettaServer.py

    Once run, a message should appear in the PyMOL terminal window indicating the server was 
    initialized successfully. 

3.  Run Rosetta application:  

    From the regular terminal, run the standalone application or other membrane framework apps
    from the command line: 

        Rosetta/main/source/bin/view_membrane_protein.linuxgccrelease -database /path/to/db @pymol_flags

    Within ~10-20 seconds, 2 parallel planes (PyMOL object entitled membrane_planes) will appear 
    in the PyMOL session. 

    Note: In earlier versions of PyMOL, small lines will be drawn to the MEM residue as this is a HETATM
    in the PDB. These can be turned off manually by hiding the line view for the MEM residue. 

Example Outputs
---------------
The following outputs will be generated from the standalone pymol viewer application. A version of these outputs 
are also provided in the example_outputs/ directory: 

1. `1c3w.pse`          : Example pymol session file including membrane planes objects
2. `1c3w.png`          : Example image of PDB 1c3w (bacteriorhodopsin) embedded in the membrane (from session file)
3. `1c3w_tr_0001.pdb`  : Output PDB file containing membrane position information in MEM residue

Note: Rosetta will also output a score file (score.sc) which is not needed for this analysis. 

Additional References
---------------------

1. Baugh EH, Lyskov S, Weitzner BD, Gray JJ (2011) Real-Time PyMOL Visualization for Rosetta and PyRosetta. PLoS ONE 6: e21931.

2. DeLano W (n.d.) The PyMOL Manual: Compiled Graphics Objects (CGOs) and Molscript Ribbons. Available: http://pymol.sourceforge.net/newman/user/toc.html.

