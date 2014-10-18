Rosetta Membrane Framework Application: Visualizing Membranes in PyMOL
===========================================================================

### About this Protocol Capture
Author: Rebecca F. Alford (rfalford12@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: October 2014

Special acknowledgement to Evan Baugh for assistance working with the PyMOL Mover. 

Rosetta Revision: #57363

Publication describing the method: 
Alford RF, Koehler Leman JK, Weitzner BD, Gray JJ (2014)
A general framework in Rosetta 3 for modeling and design of membrane proteins
PLoS ONE (in preparation)

### Description
Visualizing the geometry of the membrane is an important step in analysis of membrane models. In particular, we can answer questions about interacitons at different hydrophobic layers, membrane embedding, and orientation. This application enables visualization of output Rosetta models with the Rosetta PyMOL viewer. This tool uses the position of the bilayer described by the membrane framework to derive the position of two parallel planes separated by the membrane thickness. These planes are modeled as CGO objects in python and are drawn on the fly during simulations. This tool can also be used in real-time during simulations with any JD2-supported Rosetta application using the same set of flags. 

### Executable/Script
The PyMOL viewer with membranes can be used with any Rosetta application that uses the membrane framework with teh same flags. We also provide an application to explicitly view output models in the C++ code: 

Application: view_membrane_protein.<platform-exe> 

### Generating Inputs
Viewing the membrane planes only requires a spanfile as input. 

  1. Generating a Spanfile
  A spanfile describing transmembrane spanning regions can be generated using the OCTOPUS server
  (http://octopus.cbr.su.se/). This file must be converted to a Rosetta spanfile format using octopus2span.pl

    cd mpframework-ddG/scripts/
    ./octopus2span.pl octopus_pred.out > spanfile.txt

### Useful Scripts
This demo contains a script directory with: 
  - octopus2span.pl: Convert OCTOPUS topology prediction to Rosetta spanfile format

### Running the Application

1. Required Flags 
To run this applicaiton, the minimum required flags are described below. Flags are also included
in a flags file in this demo: 

flags                                  descriptions
--------------------------------------------------------------------------------------------------
-in:file:s <pdbfile>                   Input PDB Structure: Asymmetric input structure (should have been 
                                       generated as mystruct_input.pdb after running make_symmdef_file.pl)
-membrane_new:setup:spanfiles          Spanfile describing spanning topology of starting structure 
                                       for full symmetric structure
-show_simulation_in_pymol 0			   Use the PyMOL viewer to visualize membrane planes for structures
-keep_pymol_simulation_history 1       Keep pymol frames for making movies/replaying simulations (optional)

2. Command line

./view_membrane_protein.<exe> -database /path/to/my/rosettadb @flags

## Example Outputs
The following example outputs are included in this demo in the example_outputs/ directory: 
 - 1c3w.pse: Example pymol session file including membrane planes objects
 - 1c3w.png: Example image of PDB 1C3W (bacteriorhodopsin) embedded in the membrane (from session file)

## References
1. Baugh EH, Lyskov S, Weitzner BD, Gray JJ (2011) Real-Time PyMOL Visualization for Rosetta and PyRosetta. PLoS ONE 6: e21931.

2. DeLano W (n.d.) The PyMOL Manual: Compiled Graphics Objects (CGOs) and Molscript Ribbons. Available: http://pymol.sourceforge.net/newman/user/toc.html.
