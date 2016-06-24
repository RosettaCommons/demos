# Nucleobase Sample Around

KEYWORDS: NUCLEIC_ACIDS UTILITIES

## Author
Rhiju Das, rhiju@stanford.edu

Edited by Parisa Hosseinzadeh (parisah@uw.edu) during XRW documentation 2016 to enable automatic testing of demos.

## Brief Description

Make tables of interaction energies between an adenosine nucleobase and, say, 
 a simple carbon atom or phosphate as a probe.

## Abstract

We wanted to compare potentials of mean force of various atoms like a single methyl probe, water, adenosine, etc. around a fixed adenosine to explicit molecular dynamics solutions.

Developed in summer 2012, Das Lab hackathon -- Kyle Beauchamp, Fang-Chieh Chou, Rhiju Das, Parin Sripakdeevong. 
Extended to include phosphate by Rhiju Das in Dec. 2014.

## Running
To sample a 'carbon' probe atom:
```
 $> $ROSETTA3/bin/nucleobase_sample_around.default.linuxgccrelease  -s rosetta_inputs/a_RNA.pdb
```
where `$ROSETTA3`=path-to-Rosetta/main/source

To sample a water (sampling all possible orientations and outputting Boltzmann summed free energies)
```
 $> $ROSETTA3/bin/nucleobase_sample_around.default.linuxgccrelease  -s rosetta_inputs/a_RNA.pdb  -sample_water 
```
To sample a nucleobase
```
 $> $ROSETTA3/bin/nucleobase_sample_around.default.linuxgccrelease  -s rosetta_inputs/a_RNA.pdb  -sample_another_nucleobase   -copy_nucleobase_nucleobase_file rosetta_inputs/double_A_ready_set.pdb
```
To sample an nucleobase, reading in a starting nucleobase-nucleobase pairing conformation.
```
 $> $ROSETTA3/bin/nucleobase_sample_around.default.linuxgccrelease  -s rosetta_inputs/a_RNA.pdb  -sample_another_nucleobase   -copy_nucleobase_nucleobase_file rosetta_inputs/double_A_ready_set.pdb
```
Can now sample phosphates with the flags
```
 $> $ROSETTA3/bin/nucleobase_sample_around.default.linuxgccrelease -sample_phosphate -center_on_OP2
```
The phosphate center is on the phosphorus atom, unless user specifies -center_on_OP2 . 
Note that due to some silliness in available variant types and the desire to use a phosphate from an actual nucleotide residue_type, the probe phosphate also has a floating C1'.

Recently added option
```
 -nucleobase g [or a,c,u]
```
which will use something other than adenosine as the central nucleotide to sample around

## Plotting Results

The plotting script is available in Rosetta/tools/rna_tools/pdb_util/plot_contour.py
```
plot_contour.py score_xy_0.table score_xy_0.png
plot_contour.py score_xy_1.5.table score_xy_1.5.png
plot_contour.py score_xy_3.table score_xy_3.png
plot_contour.py score_xy_4.table score_xy_4.png
```

