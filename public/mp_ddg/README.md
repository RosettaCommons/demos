Membrane ΔΔG
============

KEYWORDS: MEMBRANES DOCKING

Author: Rebecca F. Alford (rfalford12@gmail.com)   
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)  
Last Updated: June 2016  

---

Measuring free energy changes upon mutation can inform our understanding of membrane protein stability and variation and is a step toward design. In this application, we predict ddGs by measuring the difference in Rosetta energy for the native and mutated conformation. This application uses a modified version of the all atom energy function for membrane proteins, which includes the fa_elec term and pH energy (see below). The Membrane ddG application is part of the RosettaMP Framework.

Documentation Link:  
* https://www.rosettacommons.org/docs/wiki/Membrane-ddG

Publication describing the method:  
* Alford RF, Koehler Leman J, Weitzner BD, Gray JJ (2015) An integrated 
  framework advancing membrane protein modeling and design PLosCompBio (Under 
  Review) 

## Executable/Script ##
The membrane ddG application is implemented as a python script in PyRosetta. The scripts described here can be found in this protocol capture. Developmental versions can also be found: 

1. In the Rosetta Source Code: 

        /path/to/Rosetta/source/src/python/bindings/app/membrane/predict_ddG.py

2. In the PyRosetta package: 

        /path/to/PyRosetta/app/membrane/predict_ddG.py

## Generating Inputs ##
Three inputs are required for the ddG application:  

1. PDB file for the protein structure transformed into the membrane coordinate frame.
2. Span file describing the location of trans-membrane spans
3. Residue position to mutate to (uses pose numbering)

Steps for generating these inputs are found below. A set of example inputs can 
also be found in inputs/. Here, OmpLA (PDB ID: 1qd6) is used as an example: 

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

   Once the PDB is downloaded from the PDBTM, clean the PDB and extract chain 
   C using teh clean_pdb.py script in Rosetta tools using the following command: 

        /path/to/Rosetta/tools/protein_tools/scripts/clean_pdb.py 1qd6_tr.pdb C

   The resulting PDB is referred to as 1qd6_tr_C.pdb from here. 

2. Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2015. An example commandline using 
   1qd6 is also provided here: 

        Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -in:file:s inputs/1qd6_tr_C.pdb

   For this example, this command will produce 1 output file: 
     = 1qd6_tr_C.span: Spanfile containing predicted trans-membrane spans

   Note: For this example, 1qd6 should have 12 transmembrane spans. Adjust the spanfile,
   if needed. 

## Steps for Running the protocol ##
Here, we describe the steps required to run the MP_ddG protocol. First, we describe how to 
assemble a simple PyRosetta script using the membrane framework for ddG predictions (predict_ompLA_ddG.py). Next, we describe use of a general ddG prediction application for larger scale use. 

1. Application-Specific Membrane ddG PyRosetta Protocol
   PyRosetta calculations can be adapted to use the Rosetta Membrane Framework
   with only a few additional steps. These include: 

   * Use AddMembraneMover (in protcols.membrane) to initialize the membrane framework
   * Use MembranePositionFromTopologyMover to orient the pose in the membrane based on the transmembrane spans (optional, but recommended)
   * Setup a membrane energy function
   * Proceed with normal Rosetta functions

   Here, we provide an example application-specific ddG calculation script for computing ddGs of mutation in OmpLA for comparison with experimental values in Moon & Fleming, 2011. The script can be run with no arguments by the following command: 

        ./predict_OmpLA_ddG.py 

   Step-by-step instructions on how to setup this script are provided in the predict_OmpLA_ddG.py script (in this protocol capture). 

   A single output file is created by this script: 

   * ompLA_ddG.out: Predicted ddGs for each mutation

2. Large Scale ddG predictions with the RosettaMP Framework
   Here, we describe the steps to run the MPddG protocol, incorporating both
   repacking and pH effects. Additional options are described in the documentation above. Here, we also use OmpLA as an example: 

   Here, you will need to specify the input PDB, spanfile, and residue position to 
   mutate. By default, ddGs to all canonical residues will be computed. A specific 
   ddG of mutation can be computed using the flag --mut <AA>. In this example, we specify a repack radius of 8.0A. This means all residues within 8A of the mutant position are repacked. We also specify the pH at which predictions are carried out. 

   This application can be run using the following command line: 

        ./predict_ddG.py --in_pdb inputs/1qd6_tr_C.pdb --in_span inputs/1qd6_tr_C.span --res 181 --repack_radius 8.0

   Two output files (with default names) are created by this script
   * ddG.out: predicted ddGs per mutation
   * scores.sc: Breakdown of ddGs by Rosetta score term (weighted)

## Example Outputs
Outputs from the two scripts above were renamed for clarity

1. The predict_OmpLA_ddG.py script will write a list of mutations and ddG 
   values to an output file. The columns in the file are residu position, 
   mutant amino acid (1-letter code) and predicted ddG

   An example is provided in example_outputs/OmpLA_predicted_ddGds_specific.out

2. The predict_ddG.py script will create two output files: 

   * A list of mutations and their predicted ddGs. The columns in this file are
     pdb name, mutant amino acid, mutant score, native score, ddG. An example is provided in

        example_outputs/OmpLA_predicted_ddGds_general.out

   * A breakdown of each ddG by Rosetta score term (weighted). An example 
     output file is provided in

        example_outputs/OmpLA_ddG_breakdown.sc

For all scripts - running multiple times with the same output path specified will APPEND to the file and not overwrite it

## Version Info ##
Rosetta Revision #58069 
Python Version: 2.7+  
PyRosetta Release Version: March 2015  

PyRosetta can be downloaded from http://www.pyrosetta.org. Follow the instructions
provided in the README to setup your shell environment

## Additional References ##
1. Chaudhury S, Lyskov S, Gray JJ (2010) PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta.

2. Moon CP, Fleming KG (2011) Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proc Natl Acad Sci. 

3. Kellogg, Elizabeth H., Leaver-Fay A, and Baker D. “Role of Conformational Sampling in Computing Mutation-Induced Changes in Protein Structure and Stability.” Proteins 79, no. 3 (March 2011): 830–38. doi:10.1002/prot.22921.

4. Kilambi, KP, and Gray JJ. “Rapid Calculation of Protein pKa Values Using Rosetta.” Biophysical Journal 103, no. 3 (August 8, 2012): 587–95. doi:10.1016/j.bpj.2012.06.044.

