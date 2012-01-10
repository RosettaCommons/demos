CLAIMER SequenceClaimer
LABEL main
FILE 1nsf.fasta
END_CLAIMER

CLAIMER DensityScoringClaimer
anchor 98
END_CLAIMER

CLAIMER CoordConstraintClaimer
PDB_FILE 1nsf_bad_refine.pdb
# put root into one of the helices
ASK_FOR_ROOT ALL
POTENTIAL BOUNDED 0.0 4 1 xyz
END_CLAIMER

ABINITIO_FRAGS
#this will initialize fragments mover as in classic-abinitio
#small frags and large frags   smooth moves in stage4
LARGE aansfd_09_05.200_v1_3.gz
SMALL aansfd_03_05.200_v1_3.gz
END_ABINITO

CLAIMER RigidChunkClaimer
## defines a chunk
pdb 1nsf_bad_refine.pdb
REGION
RIGID 24 42 0 0
RIGID 61 68 0 0
RIGID 90 105 0 0
RIGID 135 143 0 0
RIGID 161 171 0 0
RIGID 185 195 0 0
RIGID 200 207 0 0
RIGID 219 229 0 0
END_REGION
END_CLAIMER
