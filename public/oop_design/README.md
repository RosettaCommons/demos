# Oop Design: Design with oligooxopiperazine helix mimetics

KEYWORDS: PEPTIDES DESIGN

## Author
written by Kevin Drew (Bonneau Lab), kdrew@nyu.edu

## General Description
This demo shows how to run the oop\_design application.  An oligooxopiperazine (oop) is a helical mimetic scaffold used for inhibiting protein interactions. The demo shows the design of an oop inhibitor for the MDM2-P53 protein interaction.

## Algorithm
1. Pertubation phase: rigid body movement of oop wrt target, oop small moves to oop ring conformation, oop puck move to change oop ring pucker, small moves to ring linkers
2. Design phase: design user specified residues on oop scaffold and minimize
3. Repeat 10x

## Command
```
oop_design<.exe> @input/flags
```
for example, you can run (where `$ROSETTA3`=path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/oop_design.default.linuxgccrelease @inputs/flags
```

## Input files
```
./input/flags = user specified options
./input/mdm2_oopAAAA.pdb = input structure where target is chain 1 and oop is chain 2
```

## Options
```
-oop_design_positions = residues on oop to design (numbering is relative to oop, for example 3 is the third residue on oop), default repacks with no design
-pert_num = number of pertubations during pertubation phase, default 10, production 100
-design_loop_num = number of pertubation + design cycles, default 10, production 10
-nstruct = production 1000
```

## Post Processing
Similar to other multi chain design protocols, the ddG is computed and is a good indicator of a good design.  First sort by total score, take top 5 percent and then sort by REPACK_ENERGY_DIFF.

## Limitations: 
This app is inflexible to adjusting monte carlo temperatures, score functions, degree of rigid body pertubations, designing noncanonical amino acids, etc. The app also requires the oop is close to a plausible binding mode with respect to the target. 

