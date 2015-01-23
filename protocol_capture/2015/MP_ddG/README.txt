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

PyRosetta can be downloaded from http://www.pyrosetta.org

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

## Generating Inputs ##
Two inputs are required for the ddG application:  
  (1) PDB file for the protein structure (membrane-transformed)
  (2) Span file describing the location of trans-membrane spans

Steps for generating these inputs are found below. A set of example inputs can 
also be found in inputs/. Here, OmpLA (PDB ID: 1qd6) is used as an example: 

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2014. An example commandline using 
   1qd6 is also provided here: 

   Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/db -in:file:s inputs/1qd6_tr.pdb

   For this example, this command will produce 1 output files: 
     = 1qd6_tr.span: Spanfile containing predicted trans-membrane spans

   Note: For this example, 1qd6 should have 12 transmembrane spans

## Steps for Running each Protocol ##
Here, we describe the steps required to run the MP_ddG protocol. As an example, all steps 
use the PDB 1qd6: 

  (1) OmpLA Example: To run the membrane ddG application for the ompLA example, run the following command
      (no arguments). Inside this script is also instructions for how to include the membrane
      framework in a general PyRosetta script: 

		./compute_ompLA_ddG.py

	The ddGs are in the log output printed to the screen.

  (2) General Script: Use the general case script to run the MP_ddG application for a general case.
      Here, you will need to specify the input PDB, spanfile, output file name containing predicted ddGs, 
      residue number to mutate, and optionally the mutant you are interested in. Run the following
      command line: 

		./compute_ddG.py -p inputs/1qd6_tr.pdb -s inputs/1qd6_tr.span -o ddGs_ompLA.txt -r 181

 	  The flags are:
 	  -p <input pdb file>
	  -s <input spanfile>
	  -o <output filename with predicted ddGs>
	  -r <residue number to mutate: this is pose residue numbering, so it's always safer to renumber your input PDB starting from 1 before running this script>
	  -m <optional: one-letter-code of amino acid to mutate into. For instance: for the 'L5A' mutation it would be 'A'>

	  The columns in the output file are: amino acid, mutant score, native score, ddG. Running the script multiple times with the same output file will APPEND to the file and not overwrite it.

## Example Outputs
  1. The compute_OmpLA_ddG.py script will write a list of mutations and ddG values to standard output. 
  2. The compute_ddG.py script will write ddGs for each mutation to a file. This file, ddG_ompLA.txt, is included
     for OmpLA in example_outputs/ 

## Additional References ##
1. Chaudhury S, Lyskov S, Gray JJ (2010) PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta.

2.  Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. 

