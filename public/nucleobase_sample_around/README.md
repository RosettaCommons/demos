### General Information ##################
# Your name:
Rhiju Das
rhiju@stanford.edu

# Protocol Name:
nucleobase_sample_around

# Brief Description:

Make tables of interaction energies between an adenosine nucleobase and, say, 
 a simple carbon atom or phosphate as a probe.

# Abstract

We wanted to compare potentials of mean force of various atoms like a single methyl probe, water, adenosine, etc. around a fixed adenosine to explicit molecular dynamics solutions.

Developed in summer 2012, Das Lab hackathon -- Kyle Beauchamp, Fang-Chieh Chou, Rhiju Das, Parin Sripakdeevong. 
Extended to include phosphate by Rhiju Das in Dec. 2014.

### running #########
# Example Rosetta Command Lines:
To sample a 'carbon' probe atom:

 nucleobase_sample_around   [-s a_RNA.pdb]

To sample a water (sampling all possible orientations and outputting Boltzmann summed free energies)

 nucleobase_sample_around   [-s a_RNA.pdb]  -sample_water  [ -extra_res ~/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/water/TP3.params ]

To sample a nucleobase

 nucleobase_sample_around   [-s a_RNA.pdb]  -sample_another_nucleobase   -copy_nucleobase_nucleobase_file double_A_ready_set.pdb

To sample an nucleobase, reading in a starting nucleobase-nucleobase pairing conformation.

 nucleobase_sample_around   [-s a_RNA.pdb]  -sample_another_nucleobase   -copy_nucleobase_nucleobase_file double_A_ready_set.pdb

Can now sample phosphates with the flags

 nucleobase_sample_around -sample_phosphate [-center_on_OP2]

The phosphate center is on the phosphorus atom, unless user specifies -center_on_OP2 . 
Note that due to some silliness in available variant types and the desire to use a phosphate from an actual nucleotide residue_type, the probe phosphate also has a floating C1'.

Recently added option

 -nucleobase g [or a,c,u]

which will use something other than adenosine as the central nucleotide to sample around



### plotting results ###

The plotting script is available in Rosetta/tools/rna_tools/pdb_util/plot_contour.py

plot_contour.py score_xy_0.table score_xy_0.png
plot_contour.py score_xy_1.5.table score_xy_1.5.png
plot_contour.py score_xy_3.table score_xy_3.png
plot_contour.py score_xy_4.table score_xy_4.png
