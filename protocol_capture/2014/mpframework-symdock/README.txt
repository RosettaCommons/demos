
Rosetta Membrane Framework Application: Symmetric Protein-Protein Docking
===========================================================================

### About this Protocol Capture
Author: Rebecca F. Alford (rfalford12@gmail.com)
Corresponding PI: Jeffrey J. Gray (jgray@jhu.edu)
Last Updated: December 2014

Rosetta Revision #57514

Alford RF, Koehler Leman J, Weitzner BD, Duran A, Elazar A, Tiley D, Gray JJ (2015)
An integrated framework advancing membrane protein modeling and design
PLoS ONE (in preparation) 

### Description
Membrane proteins often assemble into symmetric complexes in the membrane environment. While the exact reason for internal symmetry is stil debated, it has been suggested that symmetry enhances proteins stability and avoid traps in the energy landscape during conformational changes needed for function. The symmetric membrane protein-protein docking application combines the symmetry machinery, membrane framework, and symmetric docking applciation to assemble symmetric complexes in the membrane. 

## Algorithm
The application first constructs a symmetric membrane protein conformation from an asymmetric unit by placing a membrane virtual at the origin of the coordinate frame and arranging subunits around the membrane normal axis (cyclic symmetry). The symmetric docking application is then used to apply symmetric moves, packing operations, and iterations of scoring using the membrane all atom energy function to sample possible conformations. The ensemble of output models describe possible conformations of this complex in the membrane environment. 

This application currently supports modeling of cyclic symmetries. 

### Executable/Script
Membrane Symmetric Protein-Protein Docking is a C++ application in the Rosetta 3 software suite. 

Application: membrane_symmdocking.<platform-exe> 

## Generating Inputs
The symmetric membrane protein-protein docking application requires two input files: 

  1. Generating a Spanfile
  A spanfile describing transmembrane spanning regions can be generated using the OCTOPUS server
  (http://octopus.cbr.su.se/). This file must be converted to a Rosetta spanfile format using the
  octopus2span.pl script. Use the following command: 

    cd mpframework-symdock/scripts/
    ./octopus2span.pl octopus_pred.out > spanfile.txt

  2. Generate a Symmetry definition file
  Symmetry definition files describe the symmetry of the system and serve as instructions to Rosetta for building the full starting complex from an asymmetric unit. This file can be generated using the 
  make_symmdef_file.pl script in src/apps/public/symmetry/. Detailed information for generation of symmetry definition files can be found at https://www.rosettacommons.org/docs/latest/make-symmdef-file.html. In general, this script requires a PDB file, symmetry mode (NCS or PSEUDO only for membranes), specification of the asymmetric unit and symmetry group # if known. 

  The symmetry definition file in this example was generated with the following commandline. 

    ./make_symmdef_file.pl -m NCS -p 1afo.pdb -r 10.0 -a A -b B:2 

  Symmetries can also be generated using a denovo protocol. If the user knows the symmetry group (e.g. C3), 
  the user can generate a definition for this symmetry group independent of the structure. This script can 
  be found in src/apps/pubilc/symmetry/make_symmdef_denovo.py. More information for using this script can 
  be found here: https://www.rosettacommons.org/docs/latest/make-symmdef-file-denovo.html

## Useful Scripts
This demo contains a script directory with: 
  - octopus2span.pl: Convert OCTOPUS topology prediction to Rosetta spanfile format
  - predict_lips.pl: Use the TMPLIP server to predict per-residue lipophilicity
  - alignblast.pl: Perform and parse multiple sequnece alignment from psiblast (needed for   
    predict_lips.pl)

The make_symmdef_file scripts can be found in src/apps/public/symmetry. 
    
### Running the Application

1. Required Flags 
To run this applicaiton, the minimum required flags are described below. Flags are also included
in a flags file in this demo: 

flags                                  descriptions
--------------------------------------------------------------------------------------------------
-in:file:s <pdbfile>                   Input PDB Structure: Asymmetric input structure (should have been 
                                       generated as mystruct_input.pdb after running make_symmdef_file.pl)
-membrane_new:setup:spanfiles          Spanfile describing spanning topology of starting structure 
                                       for full symmetric structure
-membrane_new:scoring:hbond            Turn on depth-dependent hydrogen bonding term when using the   
                                       membrane high resolution energy function
-symmetry:symmetry_definition          Symmetry definition file
-symmetry:initialize_rigid_body_dofs   Initialize rigid body configuration (symmetric)
-nstruct                               Number of structures to generate

2. Recommended # of Decoys
 - For demo run: 1
 - For production runs: 1000

3. Command line
To run this application, use the following command line: 

./membrane_symdocking.<exe> -database /path/to/my/rosettadb @flags 

### Example Outputs
The following example outputs are included in this demo in the example_outputs/ directory: 
 - 1afo_tr_out.pdb: Output decoy from symmetric docking of 1afo
 - score.sc: Scorefile output by Rosetta contianing membrane and symmetry scores 

### References 
1. DiMaio F, Leaver-Fay A, Bradley P, Baker D, André I (2011) Modeling Symmetric Macromolecular Structures in Rosetta3. PLoS ONE 6: e20450. 

2. Barth P, Schonbrun J, Baker D (2007) Toward high-resolution prediction and design of transmembrane helical protein structures. Proc Natl Acad Sci 104: 15682–15687. 
