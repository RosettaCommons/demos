MP_Framework Application: Membrane Protein-Protein Docking
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
This script docks two partners of a membrane protein inside the membrane bilayer.
The app inherits from the regular DockingProtocol and uses the membrane framework.
Before running MPdocking, both partners need to be in a single PDB file (can be
accomplished using mpdocking_setup app), transformed into the membrane, and
prepacked. To do this, use the mpdocking_prepack app - demo available. 

A note on the score function: currently this application uses its own 
scorefunction (which is a combination of Vladimirs membrane protein scorefunction
and the docking scorefunction), so PLEASE don't supply a scorefunction on the 
commandline (in fact, the app might crash if you do). Preliminary data gives 
good sampling (even below 1A RMSD to the native in many cases) for local docking and the 
scorefunction is decent but not perfect (more than half of the structures give 
funnels) with lowest energies in the 5A range. The scorefunction for this will 
be further optimized soon.

## Executable/Script
Rosetta/main/source/bin/mpdocking.macosclangrelease

## Generating Inputs 
You need:
(1) a PDB file of a native structure: only for calculating meaningful RMSDs:
    here: native.pdb

 a) Generate a PDB file where the membrane protein structure (our case 1AFO) is 
    transformed into PDB coordinates (z-axis is membrane normal). This can be done 
    either by downloading the transformed PDB directly from the PDBTM website 
    (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
    it through the PPM server (http://opm.phar.umich.edu/server.php).

 b) Clean the PDB file by using clean_pdb.pl in the folder 
    Rosetta/tools/protein_tools/scripts/:
    >>clean_pdb.pl 1AFO_tr.pdb ignorechain

(2) a PDB file transformed into membrane coordinates, cleaned, prepacked:
    here: 1AFO_AB_0001.pdb

    You can use the native from step 1. If not, make sure both partners are in a
    single PDB file, transformed into membrane coordinates, and cleaned as in step
    1. Run the MP prepacking protocol described in the demo MPdocking_prepack. 
    MAKE SURE YOU REMOVE THE MEM RESIDUES IN THE PREPACKED AND IN THE NATIVE 
    STRUCTURE!!! A future revision will have the option to read in a membrane residue
	but the current revision does not have this functionality for MP_docking. 

(3) a span file containing the topology information of the membrane protein
    here: 1AFO_AB.span

 a) Prepare the spanfile from the PDB file according to the instructions in the 
    demo spanfile_from_pdb.

Example inputs are found in the input folder. 

## Running the Application
Run the command.sh script provided in this folder:
>>./command.sh

## Example Outputs
Outputs are found in the output folder:

1) 1AFO_AB_0001_0001.pdb - final decoy with MEM residue

2) score_mpdock_1AFO.sc  - score file
   If you supplied a proper native, the rms column in the score file should be the
   RMSD to the native structure. The rms in column 3 is the ligand RMSD, i.e. the 
   RMSD over the docked partner, disregarding the fixed partner. The Irms in column
   5 is the interface RMSD.

