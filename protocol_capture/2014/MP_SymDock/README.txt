
Rosetta Membrane Framework Application: Symmetric Protein-Protein Docking
===========================================================================

### About this Protocol Capture
Author: Rebecca F. Alford (rfalford12@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: January 2015

Rosetta Revision #57518

Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ (2015)
An integrated framework advancing membrane protein modeling and design
PLoS ONE (in preparation) 

## Description ## 
This application assembles and docks symmetric protein complexes in the membrane
bilayer. The full symmetric native complex is first refined using the MP_Relax
application. The lowest scoring native conformation (by total Rosetta score) is
then used as input to the membrane symmetric docking application, which searches
for possible confromations by reassembling and docking subunits together. 

This application combines the membrane framework, symmetry machineary, and standard
symmetric docking algorithm in Rosetta. Currently, docking of Cyclic (C) symmetries is
supported. 

## Executable/Script ##
Rosetta/main/source/bin/membrane_symmdocking.linuxgccrelease

## Generating Inputs ##
Three initial input files are required for this protocol: 
  (1) PDB file for the native symmetric complex (all subunits)
  (2) Span file describing trans-membrane spans of the full complex
  (3) Span file describing trans-membrane spans of the asymmetric unit

1. PDB File: Generate a PDB file where the membrane protein structure is transformed 
   into PDB coordinates (z-axis is membrane normal). This can be done 
   either by downloading the transformed PDB directly from the PDBTM website 
   (http://pdbtm.enzim.hu/) or by downloading a PDB file from the PDB and running
   it through the PPM server (http://opm.phar.umich.edu/server.php).

2. Full & Asymmetric Span File: Generate a spanfile from the PDB structure using
   the spanfile_from_pdb application described in the MP_spanfile-from-pdb protocol
   capture in Rosetta/demos/protocol_captures/2014. An example commandline using 
   1bl8 as an example is also provided here: 

   Rosetta/main/source/bin/spanfile_from_pdb.linuxgccrelease -database /path/to/db -in:file:s 1bl8_tr.pdb

   This command will produce 5 output files: 
     = 1bl8_tr.span: Predicted trans-membrane spans for the full symmetric complex
     = 1bl8_tr<A-D>.span: Predicted trans-membrane spans for each chain in the complex

    These inputs can also be found in the example_inputs directory

## Steps of the protocol ##
Here, we describe the steps required to run the MP_SymDock protocol. All steps 
use a C4 Symmetric Potassium Channel (PDB ID: 1bl8) as an example: 

1. Initial Refinement: Using the native symmetric complex and full spanfile, generate 
   10 refined models using the MP_Relax protocol. This protocol described in the 
   a protocol capture in Rosetta/demos/protocol_capture/2014/MP_relax. Run the following
   commandline with the given flags file: 

   Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -parser:protocol membrane_relax.xml @relax_flags

   The following output files will be generated: 
     = 1bl8_tr_<1-5>.pdb    : 5 refined models of 1bl8
     = relax_scores_1bl8.sc : Rosetta scores for each resulting models

   Examples of these outputs can be found in example_refined_models

2. Input model selection: Use the score file from the refinement step to select the lowest scoring
   refined model by total Rosetta score. This structure will be used as input in the next step and is
   referred to as 1bl8_refined.pdb from this point forward.  

3. Generate inputs for symmetry: To prepare the structure for assembly and docking in the protocol, 
   a set of asymmetric inputs must be generated. These inputs describe the asymmetric unit, which will 
   later be used to re-assemble the complex based on a defined symmetry definition. A version of these
   output files are provided in the example_symmetry_files directory

   First, create the asymmetric input structure and symmetry definition file using the make_symmdef_file.pl
   script. An example commandline is provided below: 

   Rosetta/main/source/bin/make_symmdef_file.pl -p 1bl8_refined.pdb -a A -i B:4 > 1bl8.c4.symm

   In this commandline, -p specifies the input PDB file, -a specifies the chain or chains to use as the asymmetric unit, 
   and -i specifies how to organize the remaining chains. In this example, "B:4" means use chain B as the next subunit
   and arrange subunits as a C4 tetramer (4 subunits around the Z axis). 

   This command will generate various files, only two of which are key here: 
       = 1bl8_refined_INPUT.pdb  : PDB file containing the asymmetric subunit
       = 1bl8.c4.symm            : Symmetry definitoin file describing the arrangement of subunits in the complex. 
                                   This file specifically describes needed translations and rotations to regenerate 
                                   and assemble this complex from the input file. 

    Important notes: To generate a correct symmetry, make_symmdef_file.pl requires all chains be of equal length. Subunits
    shoould also be close to 0Å rmsd to one another. Any asymmetry may result in an incorrect symmetry definiton. To check, 
    you can visualize the 1bl8_refined_model.pdb to ensure this initial setup is correect. 

  Next, you will also need the 1bl8_trA.span file which contains trans-membrane spans for only the asymmetric unit. 
  If your asymmetric unit contains multiple chains, you may need to assemble this file yourself from the full set
  of spans. 

4. Running the symmetric docking application
Using the asymmetric unit PDB, symmetry definition file, and asymmetric unit span file as inputs, you are ready to run the 
membrane symmetric docking application. Flags, recommended settings, and commandlines are described below: 

  (1) Required flags: flags needed to run this application are described below. A file with these flags, symdock_flags, 
  is also provided for the 1bl8 example. 

  flags                                  descriptions
  --------------------------------------------------------------------------------------------------
  -in:file:s <pdbfile>                        Input PDB Structure: Asymmetric input structure
  -membrane_new:setup:spanfiles <spanfile>    Spanfile describing spanning topology of asymmetric unit
  -membrane_new:scoring:hbond                 Turn on depth-dependent hydrogen bonding term when using the   
                                              membrane high resolution energy function
  -symmetry:symmetry_definition               Symmetry definition file
  -symmetry:initialize_rigid_body_dofs        Locally sample rigid body conformations during intial complex assembly
                                              (before docking algorithm)
  -nstruct                                    Number of structures to generate
  -packing:pack_missing_sidechains 0          Wait to pack until the membrane mode is turned on

  (2) Recommended # of Decoys
    - For demo run: 1
    - For production runs: 1000

  (3) Command Line
  Rosetta/main/source/bin/membrane_symdocking.linuxgccrelease -database /path/to/my/rosettadb @symdock_flags 

### Example Outputs
The folowing outputs will be generated from the symmetric docking protocol. A version of these outputs are also
provided in the example_outputs/ directory: 

  1. 1bl8_tr_input_0001.pdb: Symmetrically docked output model from the protocol
  2. score.sc: Scorefile output by Rosetta containing memrbane and symmetry scores for this model

### References 
1. DiMaio F, Leaver-Fay A, Bradley P, Baker D, André I (2011) Modeling Symmetric Macromolecular Structures in Rosetta3. PLoS ONE 6: e20450. 

2. Barth P, Schonbrun J, Baker D (2007) Toward high-resolution prediction and design of transmembrane helical protein structures. Proc Natl Acad Sci 104: 15682–15687. 
