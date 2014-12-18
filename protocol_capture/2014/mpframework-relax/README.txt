Rosetta Membrane Framework Application: Membrane Relax
===========================================================================

### About this Protocol Capture
Author: Rebecca F. Alford (rfalford12@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: December 2014

Rosetta Revision #57514

Publication describing the method: 
Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ (2015)
An integrated framework advancing membrane protein modeling and design
PLoS ONE (in preparation) 

### Description
Structural refinement can reveal an ensemble of states describing the conformation of a protein. It is also necessary to advance many low resolution structures determined by x-ray crystallography to atomic-level detail. These refined structures are used as inputs to several new modeling protocols such as protein-protein docking, ligand docking, etc. 

### Algorithm Description
The membrane relax application combines the traditional fast relax algorithm with the membrane all atom energy function and an addiitonal step for optimizing the membrane position. To refine the protein structure, several iterations to sample backbone and side chain conformations. In addition, the relative orientation of the membrane and protein is sampled by minimizing the transform in the membrane jump. Combination of high-resolution refinement and orienting the membane posiiton allows for simultaneous refinement of structure and optimization of membrane embedding. 

### Executable/Script
This application uses the membrane framework and is implemented as a Rosetta script. THe script included in this file is called membrane_relax.xml which can be run via the 
rosetta_scripts executable. 

### Generating Inputs
The membrane relax application requires 1 input file: 

  1. Generating a Spanfile
  A spanfile describing transmembrane spanning regions can be generated using the OCTOPUS server (http://octopus.cbr.su.se/). This file must be converted to a Rosetta spanfile format using octopus2span.pl. Example command is given below: 

    cd mpframework-relax/scripts/
    ./octopus2span.pl octopus_pred.out > spanfile.txt

### Useful Scripts
This demo contains a script directory with: 
  - octopus2span.pl: Convert OCTOPUS topology prediction to Rosetta spanfile format
  - predict_lips.pl: Use the TMPLIP server to predict per-residue lipophilicity
  - alignblast.pl: Perform and parse multiple sequnece alignment from psiblast (needed 
    for predict_lips.pl)
    
### Running the Application

1. Required Flags 
To run this applicaiton, the minimum required flags are described below. Flags are also included
in a flags file in this demo: 

flags                                  descriptions
--------------------------------------------------------------------------------------------------
-parser:protocol membrane_relax.xml    Specify membrane relax protocol to rosetta scripts executable
-in:file:s <pdbfile>                   Input PDB Structure: Asymmetric input structure (should have been 
                                       generated as mystruct_input.pdb after running make_symmdef_file.pl)
-in:ignore_unrecognized_res            Ignore unrecognized residues during initial PDB parsing (standard
                                       Rosetta flag)
-membrane_new:setup:spanfiles          Spanfile describing spanning topology of starting structure 
-membrane_new:scoring:hbond            Turn on depth-dependent hydrogen bonding term when using the   
                                       membrane high resolution energy function
-relax:fast                            Use the FastRelax mode of Rosetta Relax (uses 5 repeat cycles)
-packing:pack_missing_sidechains false Skip initial packing (will be reapplied during membrane framework 
                                       initialization)
-nstruct                               Number of structures to generate

2. Recommended # of Decoys
 - For demo run: 1
 - For production runs: 1000

3. Command line
To run this application, use the following command line: 

./rosetta_scripts.<exe> -database /path/to/my/rosettadb @flags 

### Example Outputs
The following example outputs are included with this demo in the example_outputs/ directory: 
  -2bs2_tr_output.pdb: Output relaxed decoy
  -2bs2_score.sc: Scorefile output by relax

## Refereces
1. Tyka MD, Keedy DA, Andre I, DiMaio F, Song Y, et al. (2011) Alternate states of proteins revealed by detailed energy landscape mapping. J Mol Biol. 

2. Barth P, Schonbrun J, Baker D (2007) Toward high-resolution prediction and design of transmembrane helical protein structures. Proc Natl Acad Sci 104: 15682–15687. 

3. Fleishman SJ, Leaver-Fay A, Corn JE, Strauch E-M, Khare SD, et al. (2011) RosettaScripts: A Scripting Language Interface to the Rosetta Macromolecular Modeling Suite. PLoS ONE 6: e20161. 
