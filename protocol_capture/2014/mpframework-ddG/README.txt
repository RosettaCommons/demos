Rosetta Membrane Framework Application: Membrane ddG
===========================================================================

### About this Protocol Capture
Author: Rebecca F. Alford (rfalford12@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: December 2014

Rosetta Revision #57514

Publication describing the method: 
Alford RF, Koehler Leman J, Weitzner BD, Gray JJ (2015)
An integrated framework advancing membrane protien modeling and design
PLoS ONE (in preparation) 

### Description
Measuring free energy changes upon mutation can inform our understanding of membrane protein stability and variation and is a step toward design. In this application, we measure the difference in Rosetta energy by scoring the native and mutated proteins to compute the ddG upon mutation. This application uses the all atom energy function for membrane proteins in Rosetta with the membrane framework. 

In this protocol capture, we compute the ddG upon mutation for mutations in OmpLA described in Moon & Fleming, 2011. For each mutation, the alanine at position 181, located at the midplane of the membrane bilayer, is mutated to all 19 cannonical amino acids. The application prints each ddG as output. This PyRosetta can also be easily adapted to compute ddG for other mutations by changing the input protein and list of mutations. 

### Executable/Script
The membrane ddG application is implemented as a python script in PyRosetta. The script
included in this directory is compute_ompLA_ddG.py

## Generating Inputs
Modeling membrane proteins in Rosetta requires a spanning topology file (required). To generate a spanfile
describing the transmembrane topology, you can predict topology from seuqnece using the OCTOPUS server
(http://octopus.cbr.su.se/). This file must be converted to a Rosetta spanfile format using octopus2span.pl

    cd mpframework-ddG/scripts/
    ./octopus2span.pl octopus_pred.out > spanfile.txt

## Useful Scripts
This demo contains a script directory with: 
  - octopus2span.pl: Convert OCTOPUS topology prediction to Rosetta spanfile format

## Running the Application
To run the membrane ddG script for this example case, run the python script (no arguments)

./compute_ompLA_ddG.py

## Example Outputs
The PyRosetta script will write a list of mutations and ddG values to standard output. 

### References
1. Chaudhury S, Lyskov S, Gray JJ (2010) PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta.

2.  Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. 
