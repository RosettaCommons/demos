Membrane Relax
==============

KEYWORDS: MEMBRANES STRUCTURE_PREDICTION

Author: Rebecca F. Alford (rfalford12@gmail.com)  
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: January 2015  
Rosetta Revision #58069 

---

High-resolution refinement is key for advancing low resolution structures from x-ray
crystallography to atomic level detail. For membrane proteins, this method can also
reveal an ensemble of possible membrane embeddings: the position and orientation of 
the biomolecule with respect to the membrane bilayer. 

The membrane relax application combines the Rosetta FastRelax algorithm with the
all atom energy function for membrane proteins and a gradient-based technique 
for optimizing the membrane embedding. First, a series of small backbone moves, 
rotamer trials, and minimization are used to refine the protein structure. In addition, 
the membrane position is optimizied by minimizing the "jump" or connecting relating
the MEM residue to the biomolecule. 

Publication describing the method: 
* Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ 
  (2015) An integrated framework advancing membrane protein modeling and design 
  PLoS ONE (in preparation) 

## Executable/Script ##
The membrane framework relax application is implemented in Rosetta script. This script, 
called membrane_relax.xml is included in the main directory of this protocol capture. 

It can be run with the following executable: 

    Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease

## Generating Inputs ##
Two inputs are required for the membrane relax application: 

1. PDB for the protein structure of interest

2. Span file describing the location of trans-membrane spans

Steps for generating these inputs are found below. A set of example inputs can 
also be found in example_inputs/. Here, metarhodopsin II (PDB ID: 3pxo) is 
used as an example: 

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2014. An example commandline using 
   3pxo is also provided here: 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/db -in:file:s example_inputs/3pxo_tr.pdb

   For this example, this command will produce 1 output file: 
   * 3pxo_tr.span: Spanfile containing predicted trans-membrane spans

## Steps of the protocol ##
Here, we describe the steps required to run the MP_Relax protocol. As an example, all steps 
use the PDB 3pxo: 

1. Required Options: Options (flags) needed to run this application. A file with these flags, 
   relax_flags, is also provided for 3pxo in this demo: 

        flags                                  descriptions
        --------------------------------------------------------------------------------------------------
        -parser:protocol membrane_relax.xml    Use the membrane relax protocol Rosetta script
        -in:file:s                             Input PDB Structure: PDB file for protein structure
        -membrane_new:setup:spanfiles          Spanfile describing trans-membrane spans of the starting structure
        -membrane_new:scoring:hbond            Turn on membrane depth-dependent hydrogen bonding weight
        -relax:fast                            Use the FastRelax mode of Rosetta Relax (uses 5-8 repeat cycles)
        -relax:jump_move true                  Allow the MEM and other jumps to move during refinement
        -nstruct                               Number of structures to generate
        -packing:pack_missing_sidechains 0     Wait to pack until the membrane mode is turned on
        -out:pdb                               Output all PDB structures of refined models
        -out:file:scorefile                    Specify destination for score file

2. Recommended # of Decoys

   - For demo run: 1
   - For production runs: 1000

3. Command line: 

    To run this application, use the following command line: 

        Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -database /path/to/db @relax_flags

Note on timing: Refinement in Rosetta is a time consuming application. Depending on avaialble
computing power and size of the protein, refinement of an individual decoy can take between 10-15min
for ~200 residues and between 0.5-1.0hrs for proteins > 200 residues. 

## Example Outputs ##
The following outputs will be generated from the relax protocol. A version of these outputs are also
provided in the example_outputs/ directory: 

* 3pxo_tr_0001.pdb      : Output refined model of 3pxo
* relax_scores_3pxo.sc  : Rosetta scores (including membrane scores) for refinement run

## References
1. Tyka MD, Keedy DA, Andre I, DiMaio F, Song Y, et al. (2011) Alternate states of proteins revealed by detailed energy landscape mapping. J Mol Biol. 

2. Barth P, Schonbrun J, Baker D (2007) Toward high-resolution prediction and design of transmembrane helical protein structures. Proc Natl Acad Sci 104: 15682–15687. 

3. Fleishman SJ, Leaver-Fay A, Corn JE, Strauch E-M, Khare SD, et al. (2011) RosettaScripts: A Scripting Language Interface to the Rosetta Macromolecular Modeling Suite. PLoS ONE 6: e20161. 
