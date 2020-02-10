# StepWise Monte Carlo of an RNA loop

KEYWORDS: LOOPS NUCLEIC_ACIDS RNA DESIGN 

## Author
Rhiju Das, rhiju@stanford.edu

Edited in 2016 by Parisa Hosseinzadeh to enable automatic testing of demos.

## Brief Description

Solve structure of an RNA loop or motif in the context of a starting structure.

## Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


## Running

### Example Rosetta Command Line (`$ROSETTA3`= path-to-Rosetta/main/source)
```
$> $ROSETTA3/bin/stepwise.default.linuxgccrelease -in:file:fasta rosetta_inputs/1zih.fasta -s rosetta_inputs/start_helix.pdb  -out:file:silent swm_rebuild.out -extra_min_res 4 9 -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior
```
Currently, we are mainly using a scorefunction with a more stringent torsional and repulsive potential, enabled by flags `-score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj`. 

To get out models (here from a test case). You can see one example of structure in 1zih_RNA.pdb in the rosetta_inputs directory.

```
$> $ROSETTA3/bin/extract_pdbs.default.linuxgccrelease -silent rosetta_inputs/swm_rebuild.out 
```

(Or use extract_lowscore_decoys.py which can be installed via tools/rna_tools/.)

### Example Rosetta Command Line for Design
Simply use a fasta file that has n's at positions you want to design.

```
$> $ROSETTA3/bin/stepwise.default.linuxgccrelease -in:file:fasta rosetta_inputs/NNNN.fasta -s rosetta_inputs/start_helix.pdb  -out:file:silent swm_design.out -extra_min_res 4 9 -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior
```

