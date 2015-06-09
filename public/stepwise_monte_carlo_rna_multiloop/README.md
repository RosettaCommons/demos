### General Information ##################
# Your name:
Rhiju Das
rhiju@stanford.edu

# Protocol Name:
StepWise Monte Carlo (examples for RNA)

# Brief Description:

Solve structure of an RNA internal loop or multi-helix junction.

# Abstract

Ab initio and comparative modeling of biopolymers (RNA, protein, protein/RNA) often involves solving well-defined small puzzles (4 to 20 residues), like RNA aptamers, RNA tertiary contacts, and RNA/protein interactions. If these problems have torsional combinations that have not been seen previously or are not captured by coarse-grained potentials, most Rosetta approaches will fail to recover their structures.  This app implements a stepwise ansatz, originally developed as a 'stepwise assembly' enumeration that was not reliant on fragments or coarse-grained modeling stages, but was computationally expensive. The new mode is a stepwise monte carlo, a stochastic version of stepwise assembly. 


### running #########
# Example Rosetta Command Line:

Following is for an internal loop ('two-way junction') drawn from the most conserved domain of the signal recognition particle, a core component of the machinery that translates membrane proteins in all kingdoms of life.

If you do not know the rigid body orientations of two helices (typical use case), run:

stepwise -in:file:fasta rosetta_inputs/1lnt.fasta -s rosetta_inputs/gu_gc_helix.pdb  rosetta_inputs/uc_ga_helix.pdb -out:file:silent swm_rebuild.out -extra_min_res 2 15 7 10 -terminal_res 1 8 9 16 -nstruct 20  -cycles 1000  -score:rna_torsion_potential RNA11_based_new  -native rosetta_inputs/native_1lnt_RNA.pdb

If you have starting coordinates for the two helix endpoints, you can start with that single PDB ('start_native_1lnt_RNA.pdb') instead:

stepwise -in:file:fasta rosetta_inputs/1lnt.fasta -s rosetta_inputs/start_native_1lnt_RNA.pdb -out:file:silent swm_rebuild.out -extra_min_res 2 15 7 10 -terminal_res 1 8 9 16 -nstruct 20  -cycles 1000  -score:rna_torsion_potential RNA11_based_new  -native rosetta_inputs/native_1lnt_RNA.pdb

To get out models:

extract_pdbs -silent swm_rebuild.out 

(Or use extract_lowscore_decoys.py which can be installed via tools/rna_tools/.)





