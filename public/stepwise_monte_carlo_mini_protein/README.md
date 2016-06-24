# StepWise Monte Carlo to build a small protein

KEYWORDS: STRUCTURE_PREDICTION DENOVO

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Solve structure of a mini-protein

## Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running
### Example Rosetta Command Line

```
$> $ROSETTA3/bin/stepwise.default.linuxgccrelease -fasta rosetta_inputs/2jof.fasta -native rosetta_inputs/2jof.pdb -score:weights stepwise/protein/protein_res_level_energy.wts -silent swm_rebuild.out -cycles 2 -nstruct 5
```

Most of the simulation may be spent flickering bits of secondary structure -- in the future, we will probably setup some precomputation of these bits so that computation can be focused on build up of the complete mini-protein structure.



