# Stepwise Monte Carlo RNA Multiloop

KEYWORDS: NUCLEIC_ACIDS LOOPS RNA STRUCTURE_PREDICTION

## Authors
Rhiju Das, rhiju@stanford.edu

Updated in 2016 by Parisa Hosseinzadeh (parisah@uw.edu) to enable automatic testing of demos.

## StepWise Monte Carlo (examples for RNA)

### Brief Description

Solve structure of an RNA internal loop or multi-helix junction.

### Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running

Following is for an internal loop ('two-way junction') drawn from the most conserved domain of the signal recognition particle, a core component of the machinery that translates membrane proteins in all kingdoms of life.

If you do not know the rigid body orientations of two helices (typical use case), run: (`$ROSETTA3`= path-to-Rosetta/main/source)

```
$> $ROSETTA3/bin/stepwise.default.linuxgccrelease @uk_orientation.options -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior
```

If you have starting coordinates for the two helix endpoints, you can start with that single PDB ('start_native_1lnt_RNA.pdb') instead:

```
$> $ROSETTA3/bin/stepwise.default.linuxgccrelease @known_ends.options -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior
```
For the purpose of demo, we have lowered the number of generated structures and the cycles, but usually you want to at least run 1000 cycles and generate more structures.

To get out models (in this case from the pre-generated file in the rosetta_inputs directory):

```
$> $ROSETTA3/bin/extract_pdbs.default.linuxgccrelease -silent rosetta_inputs/swm_rebuild.out -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior
```

(Or use extract_lowscore_decoys.py which can be installed via tools/rna_tools/.)

