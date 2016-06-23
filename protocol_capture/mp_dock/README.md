Membrane Protein-Protein Docking
================================

KEYWORDS: MEMBRANES DOCKING

Author: Julia Koehler Leman (julia dot koehler1982 at gmail dot com)  
Corresponding PI: Jeffrey J. Gray (jgray at jhu dot edu)  
Last Updated: January 2015  
Rosetta Revision #58069

---

To run the full application, both steps for prepacking and actual docking needs
to be carried out. 

Prepacking Step
---------------

This script prepacks a membrane protein structure inside the membrane to prepare
it for membrane protein docking. It pull the partners apart along the axis between
the center-of-masses projected onto the membrane plane, does a side-chain packing
round, and pulls them back together.

#### Executable/Script

    Rosetta/main/source/bin/docking_prepack_protocol.macosclangrelease

#### Generating Inputs

You need:

1. A PDB file transformed into membrane coordinates and cleaned

2. A span file containing the topology information of the membrane protein

PDB:

1. Generate a PDB file where the membrane protein structure (our case 1AFO) is 
   transformed into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Clean the PDB file by using clean_pdb.pl in the folder 
   Rosetta/tools/protein_tools/scripts/:
   $ clean_pdb.pl 1AFO_tr.pdb ignorechain

Span file:

1. Prepare the spanfile from the PDB file according to the instructions in the 
   demo spanfile_from_pdb

#### Running the Application

Adjust the paths in the command.sh script provided in this folder and run it:

    $ ./command.sh

For the paper, we used 10 output models and selected the lowest scoring model by
total Rosetta score. More models might be better, even though probably not totally
necessary since this does mainly sidechain optimization.

#### Example Outputs

Outputs are found in the output folder:

1. 1AFO_AB_0001.pdb  - final decoy

   The model should have MEM residues at the end (HETATM lines at residue number 81) 
   representing the virtual membrane residue. This needs to be removed for mpdocking in
   this revision even though for newer revisions you can leave that in.

2. score_ppk_1AFO.sc - score file

   Total scores can be >0 because of clashes. This is not worrysome because you 
   should run membrane protein docking afterwards which should get rid of the clashes.

MP Dock
-------

This script docks two partners of a membrane protein inside the membrane bilayer.
The app inherits from the regular DockingProtocol and uses the membrane framework.
Before running MPdocking, both partners need to be in a single PDB file (can be
accomplished using mpdocking_setup app), transformed into the membrane, and
prepacked. To do this, use the mpdocking_prepack app - demo available. 

A note on the score function: currently this application uses its own 
scorefunction (which is a combination of the membrane protein scorefunction 
[low-resolution from Vladimir Yarov-Yarovoy and high-resolution from Patrick Barth]
and the docking scorefunction), so PLEASE don't supply a scorefunction on the 
commandline (in fact, the app might crash if you do). For updates on this, please 
check the application documentation at

https://www.rosettacommons.org/docs/latest/Application-Documentation.html

Preliminary data gives good sampling (even below 1A RMSD to the native in many 
cases) for LOCAL docking and the scorefunction is decent but not perfect 
(more than half of the structures give funnels) with lowest energies in 
the 5A range. The scorefunction for this will be further optimized soon.

#### Executable/Script

    Rosetta/main/source/bin/mpdocking.macosclangrelease

#### Generating Inputs 

You need:

1.  A PDB file of a native structure: only for calculating meaningful RMSDs:
    here: native.pdb

    * Generate a PDB file where the membrane protein structure (our case 1AFO) 
      is transformed into PDB coordinates (z-axis is membrane normal). This can 
      be done either by downloading the transformed PDB directly from the PDBTM 
      website (http://pdbtm.enzim.hu/) or by downloading a PDB file from the 
      PDB and running it through the PPM server 
      (http://opm.phar.umich.edu/server.php).

    * Clean the PDB file by using clean_pdb.pl in the folder 
      Rosetta/tools/protein_tools/scripts/:

        $ clean_pdb.pl 1AFO_tr.pdb ignorechain

2.  A PDB file transformed into membrane coordinates, cleaned, prepacked:
    here: 1AFO_AB_0001.pdb

    You can use the native from step 1. If not, make sure both partners are in a
    single PDB file, transformed into membrane coordinates, and cleaned as in step 1.
    Run the MP prepacking protocol described in the demo MPdocking_prepack. 
    If things crash (they really shouldn't), remove the MEM residue in the 
    prepacked and in the native structure. 

3.  A span file containing the topology information of the membrane protein
    here: 1AFO_AB.span

    * Prepare the spanfile from the PDB file according to the instructions in 
      the demo spanfile_from_pdb.

Example inputs are found in the input folder.

#### Running the Application

Adjust the paths in the command.sh script provided in this folder and run it:

    $ ./command.sh

For production runs, build at least 1000 models. 

#### Example Outputs

Outputs are found in the output folder:

1. 1AFO_AB_0001_0001.pdb - final decoy with MEM residue

2. score_mpdock_1AFO.sc  - score file

   If you supplied a proper native, the rms column in the score file should be the
   RMSD to the native structure. The rms in column 3 is the ligand RMSD, i.e. the 
   RMSD over the docked partner, disregarding the fixed partner. The Irms in column
   5 is the interface RMSD.

References
----------

* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tilley D, Gray JJ 
  (2015): An integrated framework advancing membrane protein modeling and 
  design, PLosONE (in preparation)
