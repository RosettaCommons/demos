# C++ RosettaSurface

KEYWORDS: SURFACES DOCKING

## Authors
Written by Michael Pacella (Graylab), mpacella88@gmail.com.  Edited by Vikram K. Mulligan (vmullig@uw.edu) as part of the 2016 Documentation XRW.


## General Description
This demo will describe how to run the C++ version of the RosettaSurface 
algorithm.  The test case shown here will analyze the adsorption of the 
model peptide LK-alpha to the 104 surface of calcite

## Algorithm
Simultaneous optimization of protein rigid-body orientation, backbone and 
side chain conformations on a solid surface.  

## Commands

```bash
$> $ROSETTA3/bin/surface_docking.default.linuxgccrelease @rosetta_inputs/flags
```

In the above, ".default.linuxgccrelease" may need to be updated for your build, operating system, and compiler.

## Input Files
- lk_alpha_calcite.pdb = input pdb with LK alpha positioned above calcite 104
- calcite.surf = file containing surface vectors specific to the calcite surface in the input pdb
- lk_alpha_3mers<9mers> = 3mer and 9mer fragment files for LK alpha
- flags = arguments for rosetta surface

## Pre-processing
1. Ensure that the input pdb is properly formatted with the protein appearing before the surface in 
the input file and belonging to a separate chain

2.  Ensure that the specified surface vectors match the input surface

3.  Ensure that parameter files corresponding to the molecules comprising the surface exist in 
the rosetta database

## Post-processing

1. GetTop.sh Ads 4
2. cd TOP4.Ads
3. PostProcessRS.sh
    
Make sure that the folder contains either only adsorbed state PDBs (and native.pdb) or solution state PDBs. 
Also make sure that all post processing scripts (found in the scripts/post_processing directory)
are present in this directory as well 

These commands will extract the top 4 adsorbed-state decoys for analysis and generate
secondary structure, protein-protein contact maps, and protein-surface contact maps

## Example output
- Ads_SecStruct.png = adsorbed state secondary structure
- Sol_SecStruct.png = solution state secondary struccture
- Ads_ContactMap = adsorbed state contact map
- Sol_ContactMap = solution state contact map
- Surface_ContactMap.png = protein-surface contact map
- Diff_ContactMap.png = difference map between contacts in the ads/sol state
- lk_alpha_docked_to_calcite_0001.pdb = output final decoy
- score.sc = scorefile
- SolState_lk_alpha_docked_to_calcite_0001.pdb = solution state pdb
- Surface_lk_alpha_docked_to_calcite_0001.pdb = adsorbed state pdb

## Limitations
This app requires a single protein positioned above a solid surface whose parameters are 
present in Rosetta.  The protein needs to appear before the surface in the pdb file
and the two need to be separate chains

