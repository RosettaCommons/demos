MP_Framework Application: Prepack for Membrane Protein-Protein Docking
===========================================================================

## About this Protocol Capture
Author: Julia Koehler Leman (julia.koehler1982@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: 12/18/2014

Rosetta Revision #57518

Reference:
Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tilley D, Gray JJ (2015):
An integrated framework advancing membrane protein modeling and design, PLosONE (in preparation)

Documentation Link: <@jkleman insert the link here>  

## Description
This script prepacks a membrane protein structure inside the membrane to prepare
it for membrane protein docking. It pull the partners apart along the axis between
the center-of-masses projected onto the membrane plane, does a side-chain packing
round, and pulls them back together.

## Executable/Script
Rosetta/main/source/bin/docking_prepack_protocol.macosclangrelease

## Generating Inputs 
You need:
(1) a PDB file transformed into membrane coordinates and cleaned
(2) a span file containing the topology information of the membrane protein

PDB:
1) Generate a PDB file where the membrane protein structure (our case 1AFO) is 
   transformed into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2) Clean the PDB file by using clean_pdb.pl in the folder 
   Rosetta/tools/protein_tools/scripts/:
   >>clean_pdb.pl 1AFO_tr.pdb ignorechain

Span file:
1) Prepare the spanfile from the PDB file according to the instructions in the 
   demo spanfile_from_pdb

## Running the Application
Run the command.sh script provided in this folder:
>>./command.sh

## Example Outputs
Outputs are found in the output folder:

1) 1AFO_AB_0001.pdb  - final decoy 
   The model should have MEM residues at the end representing the virtual membrane
   residue. This needs to be removed for mpdocking!!!

2) score_ppk_1AFO.sc - score file
   Total scores can be >0 because of clashes. This is not worrysome because you 
   should run membrane protein docking afterwards which should get rid of the clashes.

