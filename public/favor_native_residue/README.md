#Favor Native Residue

KEYWORDS: DESIGN UTILITIES

This document briefly demos how to use favor_native_residue for Rosetta 3.4 in 
a fixed backbone context.  favor_native_residue is a mode (which debuted in 
Rosetta 2) that pushes the packer to prefer NOT making a mutation when the 
mutation is only weakly beneficial.  This is useful in design contexts where 
few mutations are desired.

The fixbb executable does not support favor_native_residue because the 
implementation is a little tricky.  (It could, but it currently does not).  The 
problem is that favor_native_residue uses constraints to do the energy 
favoring.  Constraints must be added to the Pose and ScoreFunction.  JD2 
doesn't instantiate the Pose until after the outer setup of fixbb has already 
completed, so there is no chance to load the constraints into fixbb.

To run the demo:

    $> $ROSETTA3/bin/rosetta_scripts.default.linuxgccrelease @inputs/options 
    
where `$ROSETTA3`=path-to-Rosetta/main/source

The options file is annotated.

The XML file (inputs/favor_native_residue.xml) contains two movers.  
FavorNativeResidue, applied before the packing mover, causes the 
favor-native-residue behavior.  "bonus" is the energy bonus to assign to native 
residues.  1.5 will likely be overwhelming; 0.5 will be useful; 0.05 will be 
very weak.  You'll have to tune it to your application.
