# Symmetry Examples

KEYWORDS: CORE_CONCEPTS SYMMETRY STRUCTURE_PREDICTION DESIGN DOCKING COMPARATIVE_MODELING  

## Authors
Frank Dimaio and Ingemar André
dimaio@u.washington.edu and ingemar.andre@biochemistry.lu.se

Update: 2016-22-06 
Sebastian Rämisch
raemisch@scripps.edu

## Brief Description
Examples of how to run symmetry-enabled rosetta protocols.

## Abstract

Symmetric protein assemblies play important roles in many biochemical processes. However, the large size of such systems is challenging for traditional structure modeling methods. This paper describes the implementation of a general framework for modeling arbitrary complex symmetries in Rosetta3.  We describe the various types of symmetries relevant to the study of protein structure that may be modeled using Rosetta’s symmetric framework.  We then describe how this symmetric framework is efficiently implemented within Rosetta, which restricts the conformational search space by sampling only symmetric degrees of freedom, and explicitly simulates only a subset of the interacting monomers.  Finally, we describe structure prediction and design applications that utilize the Rosetta3 symmetric modeling capabilities, and provide a guide to running simulations on symmetric systems.

## Software

Rosetta can be downloaded at <http://www.rosettacommons.org/software>

## Documentation

The applications that are exemplified are fully documented in the regular Rosetta documentation. For a general description see:  
[Symmetry User Guide](https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/symmetry)

[Documentation for symmetric docking](https://www.rosettacommons.org/docs/latest/application_documentation/docking/sym-dock)
[Documentation for fold-and-dock](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/fold-and-dock)
[Documentation for fixed backbone design](https://www.rosettacommons.org/docs/latest/application_documentation/design/fixbb)
[Documentation for comparative modeling](https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/RosettaCM)
[Documentation for symmetry definitions](https://www.rosettacommons.org/docs/latest/application_documentation/utilities/make-symmdef-file-denovo)

## Examples

There is a more detailed tutorial in `demos/public`, called **symmetric_docking_insulin_trimer_of_dimers**  
Defininately, take a look at that to understand how to do symmeric modeling in Rosetta. 
###1. Comparative modeling

(where `ROSETTA3`=path-to-Rosetta/main/source)

```
    $> $ROSETTA3/bin/minirosetta.default.linuxgccrelease @comparative_modeling/input_files/comparative_modeling-flags
```

###2. Fixed Backbone Design

```
$> $ROSETTA3/bin/fixbb.default.linuxgccrelease @fixbb/input_files/fixbb-symmetry-flags
```

###3. Fold-and-Dock
```
$> $ROSETTA3/bin/minirosetta.default.linuxgccrelease @fold-and-dock/input_files/fold-and-dock-flags
```

###4. Symmetric Docking
```
$> $ROSETTA3/bin/SymDock.default.linuxgccrelease  @symmetric_docking/input_files/symmetric-docking-flags
```
