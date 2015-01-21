Rosetta Membrane Framework Application: Membrane ddG
===========================================================================

### About this Protocol Capture
Author: Rebecca F. Alford (rfalford12@gmail.com)
Author: Julia Koehler Leman (julia.koehler1982@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: January 2015

Rosetta Revision #57514
Python Version 2.7+
PyRosetta Release Version: January 2015

Publication describing the method: 
Alford RF, Koehler Leman J, Weitzner BD, Gray JJ (2015)
An integrated framework advancing membrane protein modeling and design
PLoS ONE (in preparation) 

## Description ##
Measuring free energy changes upon mutation can inform our understanding of membrane
protein stability and variation. It is also an important step toward predicting
novel membrane protein folds and rational design. 

In this application, we measure the difference in Rosetta energy by scoring the native 
and mutated proteins to compute the ddG upon mutation. This application uses the all atom 
energy function for membrane proteins in Rosetta with the membrane framework. 

## Executable/Script ##
The membrane ddG application is implemented as a python script in PyRosetta. We provide 
two scripts: 

  1. compute_ompLA_ddG.py : This script demonstrates how to assemble a simple PyRosetta 
     for specialized modeling tasks with the membrane framework
  2. compute_ddG.py : A general version of the ddG application that will compute 
     a ddG given any inputs

## Generating Inputs
Modeling membrane proteins in Rosetta requires a spanning topology file (required). To generate a spanfile
describing the transmembrane topology, you can predict topology from seuqnece using the OCTOPUS server
(http://octopus.cbr.su.se/). This file must be converted to a Rosetta spanfile format using octopus2span.pl

    cd MP_ddG/scripts/
    ./octopus2span.pl octopus_pred.out > spanfile.txt

## Useful Scripts
This demo contains a script directory with: 
  - octopus2span.pl: Convert OCTOPUS topology prediction to Rosetta spanfile format

## Running the Application
To run the membrane ddG script for this example case, run the python script (no arguments)

./compute_ompLA_ddG.py

The ddGs are in the log output printed to the screen.
Alternatively, for a more general case, you can run the script

./compute_ddG.py

where you can specify the input PDB, input spanfile, output file name with predicted ddGs, the residue number to mutate, and optionally the mutant that you are interested in. Use the flag -h for input options. An example is the commandline:

./compute_ddG.py -p inputs/1qd6_tr.pdb -s inputs/1qd6_tr.span -o ddGs_ompLA.txt -r 181

The flags are:
-p <input pdb file>
-s <input spanfile>
-o <output filename with predicted ddGs>
-r <residue number to mutate: this is pose residue numbering, so it's always safer to renumber your input PDB starting from 1 before running this script>
-m <optional: one-letter-code of amino acid to mutate into. For instance: for the 'L5A' mutation it would be 'A'>

The columns in the output file are: amino acid, mutant score, native score, ddG. Running the script multiple times with the same output file will APPEND to the file and not overwrite it.

## Example Outputs
The PyRosetta script will write a list of mutations and ddG values to standard output. 

### References
1. Chaudhury S, Lyskov S, Gray JJ (2010) PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta.

2.  Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. 

