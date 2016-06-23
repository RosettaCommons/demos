# StepWise Monte Carlo (example for protein loop)

KEYWORDS: STRUCTURE_PREDICTION LOOPS

## Author
Rhiju Das, rhiju@stanford.edu

## Brief Description

Solve structure of a protein loop

## Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running

### Example Rosetta Command Line

```
$> $ROSETTA3/stepwise.linuxgccrelease -s rosetta_inputs/noloop_mini_1alc_H.pdb -fasta rosetta_inputs/mini_1alc.fasta -native rosetta_inputs/mini_1alc.pdb -score:weights stepwise/protein/protein_res_level_energy.wts -silent swm_rebuild.out -from_scratch_frequency 0.0 -allow_split_off false -cycles 2 -nstruct 2
```

Most of the simulation may be spent flickering bits of secondary structure -- in the future, we will probably setup some precomputation of these bits so that computation can be focused on build up of the complete mini-protein structure.



