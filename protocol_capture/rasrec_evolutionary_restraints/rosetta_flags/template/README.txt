README.txt

To start a Rasrec run using these files, you need to change the following lines:

flags_denovo:
FRAGS_3MERS - frag file containing 3mers
FRAGS_9MERS - frag file containing 9mers
FASTA_SEQUENCE - fasta sequence file

flags_rasrec:
REFERENCE_STRUCTURE - PDB File of reference structure
REFERENCE_RESIDUES - optional: rigid file containing residues for RMSD calculation
BROKER_FILE - either setup_init.tpb (Core Run) or setup_rerun.tpb (Refinement Run)

setup_init.tpb:
RESTRAINTS_FIRSTRUN - cst file containint distance restraints

setup_rerun.tpb:
CONVERGED_DISTANCES - *_converged_distances.cst
FILTERED_CONTACTMAPS - *.filtered_contactmaps.cst

The runtime of RASREC can be adjusted by changing the poolsize in flags_iterative.