#USE INPUT POSE

CLAIMER SequenceClaimer
CMD_FLAG 
END_CLAIMER

CLAIMER MembraneTopologyClaimer
END_CLAIMER

CLAIMER RigidChunkClaimer
NO_USE_INPUT_POSE
PDB ./input_files/2Z73A.pdb
REGION_FILE ./input_files/2Z73A_core.rigid
END_CLAIMER
