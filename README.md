Introduction
============
These demos are designed to guide users through sample procedures in computational modeling from the point of view of solving a specific problem. Each one is self-contained, containing all of the input files and scripts necessary to set up a rosetta run or runs -- inputs are in the starting_inputs/, scripts in scripts/. It also contains any inputs which one would prepare along the way before an actual rosetta run -- for example if in input pdb must be edited before use in rosetta -- these are in rosetta_inputs/ directory.  
Each directory contains a README file with an introduction to the task at hand and detailed step-by-step instructions.

The public/ directory contains demos that are in final form ready for use.
The pilot/ directory contains provisional or not-yet-completed demos that may not be ready for distribution.

Brief Descriptions
------------------
###*abinitio*###
ab initio structure prediction from a sequence with already-prepared fragements and secondary structure prediction.

###*AnchoredDesign*###
Design of a protein-protein interface around a known interaction (the anchor) at the same interface of one partner.

###*beta_strand_homodimer_design*###
Identify exposed (unpaired) beta-strands and design a homodimer using the unpaired strand.

###*darc*###
DARC is designed to dock small molecules in protein pockets.

###*electron_density*###
1. ***cryoEM tutorial*** -- use Rosetta to model structures into cryoEM data.
2. ***molecular replacement*** -- use Rosetta to solve difficult molecular replacement problems.

###*FlexPepDock_Refinement*###
Flexible peptide docking starting from an unrefined model.

###*fragment_picker*###
Prepare fragements for loop modeling or ab initio modeling.

###*ligand_dock*###
Docking a ligand into a protein. Separate ligand_setup and ligand_dock demos.

###*peptide_specificity*###
Peptide specificity / flexible backbone design. Designs a peptide to bind around an "anchor" residue on the target protein.

###*python*###
python implementation of relax protocol, to minimize the energy of a structure.

###*rna*###
How to prepare an RNA for Rosetta.

###*RNA_Assembly*###
Create a model of a large RNA from secondary structure prediction and limited constraint data (e.g. from mutations)

###*RNA_Denovo*###
Predict RNA 3 dimensional structure from sequence.

###*RNA_Design*###
Design minimum-energy RNA sequence for a given structure.

###*rosetta*###
Fragments, pdb, refold and score.

###*rosetta_scripts*###
A number of useful xml scripts

###*serialization*###

###*utility*###

###*zipper.sh*###
A file compression script
