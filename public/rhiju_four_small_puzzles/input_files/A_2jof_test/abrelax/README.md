KEYWORDS: STRUCTURE_PREDICTION GENERAL
AbinitioRelax.macosgccrelease -database ~/minirosetta_database  -fasta 2jof.fasta -native 2jof.pdb -frag3  aat000_03_05.200_v1_3.txt  -frag9 aat000_09_05.200_v1_3.txt -out:file:silent 2jof_abrelax.out -out:file:silent_struct_type binary -abinitio:relax -nstruct 1 -ex1 -ex2 -extrachi_cutoff 0

relax.macosgccrelease -s idealize_2jof.pdb -out:file:silent 2jof_nativerelax.out -out:file:silent_struct_type binary -database ~/minirosetta_database  -frag3  aat000_03_05.200_v1_3.txt  -frag9 aat000_09_05.200_v1_3.txt -native 2jof.pdb -nstruct 1  -ex1 -ex2 -extrachi_cutoff 0


