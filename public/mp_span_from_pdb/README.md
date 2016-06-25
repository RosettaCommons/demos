Spanfile from PDB
=================

KEYWORDS: MEMBRANES ANALYSIS

Author: Julia Koehler Leman (julia dot koehler1982 at gmail dot com)  
Corresponding PI: Jeffrey J. Gray (jgray at jhu dot edu)  
Last Updated: 12/18/2014  
Rosetta Revision #58069

---

This script generates a spanfile from a PDB file. A span file is a topology file
read into Rosetta and is generally generated from a membrane proteins sequence 
using the OCTOPUS server (http://octopus.cbr.su.se/) which is then converted into
a spanfile using octopus2span.pl. This app generates a spanfile from a PDB file
for easy testing of membrane applications, when the structure is known.

Reference:
* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tilley D, Gray JJ 
  (2015): An integrated framework advancing membrane protein modeling and 
  design, PLosONE (in preparation)

## Executable/Script

    main/source/bin/spanfile_from_pdb.macosclangrelease

## Generating Inputs 

1. Generate a PDB file where the membrane protein structure (our case 1AFO) is 
   transformed into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php). If you don't have it, you can copy the provided file:

        $> cp input/1AFO_tr.pdb .

2. Clean the PDB file by using clean_pdb.py in the folder 
   Rosetta/tools/protein_tools/scripts/:

        $> $ROSETTA_TOOLS/protein_tools/scripts/clean_pdb.py 1AFO_tr.pdb ignorechain

3. An example input is provided in the input folder: 1AFO_AB.pdb (your file will be named 1AFO_tr_ignorechain.pdb that you can rename).

## Running the Application
Run the command.sh script provided in this folder: (where `$ROSETTA3`=path-to-Rosetta/main/source)

    $> $ROSETTA3/bin/spanfile_from_pdb.default.linuxclangrelease -in:file:s input/1AFO_AB.pdb

## Example Outputs
This generates three spanfiles in the input (! unfortunately) folder. For this 
demo, these files have been moved into the output folder. 

    1AFO_AB.span   spanfile of full PDB
    1AFO_ABA.span  spanfile of chain A
    1AFO_ABB.span  spanfile of chain B in isolation (i.e. residue numbering starts from 1)

Please check the span file for errors and report errors to 
julia dot koehler1982 at gmail dot com!
